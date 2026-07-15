import os
from typing import List, Dict, Optional, Sequence, Tuple, Any
from concurrent.futures import ThreadPoolExecutor, as_completed
import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import glob
import json
from loguru import logger as logger
# logger.remove()
# logger.add(sys.stderr, level="DEBUG")

# Resolve anndata's backed sparse-dataset wrapper across versions. It lets us
# read only the requested rows of an on-disk layer instead of materializing the
# whole layer into memory (which `AnnData.layers[name]` does in backed mode).
try:  # anndata >= 0.11
    from anndata.io import sparse_dataset as _sparse_dataset
except Exception:  # pragma: no cover - older anndata
    try:  # anndata 0.8 - 0.10
        from anndata.experimental import sparse_dataset as _sparse_dataset
    except Exception:  # pragma: no cover
        _sparse_dataset = None


def _matrix_row_source(a: anndata.AnnData, layer: Optional[str]):
    """Return an on-disk, row-sliceable matrix for a backed `AnnData`.

    For ``layer is None`` this is the backed ``X``. For a named layer, the
    underlying HDF5 element is wrapped so that ``src[rows, :]`` reads only the
    requested rows from disk; accessing ``a.layers[name]`` directly would load
    the entire layer into memory. Returns a backed sparse dataset, an h5py
    dataset (dense), or — as a safe fallback for uncommon encodings — an
    in-memory matrix. All support ``src[sorted_rows, :]``.
    """
    if layer is None:
        return a.X
    elem = a.file["layers"][layer] #a.file is an h5py.File object
    if isinstance(elem, h5py.Group):
        enc = elem.attrs.get("encoding-type", "")
        if _sparse_dataset is not None and enc == "csr_matrix":
            return _sparse_dataset(elem)
        # csc / unknown sparse encoding: fall back to a full in-memory read
        # (still correct, just not row-streamed).
        return a.layers[layer]
    # Dense layer stored as an h5py dataset: fancy row indexing works directly.
    return elem


def _close_backed(a: anndata.AnnData) -> None:
    """Close a backed `AnnData`'s underlying HDF5 handle, across versions."""
    for attr in ("file", "_file"):
        try:
            getattr(a, attr).close()
            return
        except Exception:
            continue


def _read_rows_in_order(mat, row_key: np.ndarray):
    """Read ``mat[row_key, :]`` using sorted disk access, restoring row order.

    Sorting the positions makes the HDF5 reads sequential; the inverse
    permutation then restores the caller's requested order. Returns whatever
    type ``mat`` yields (sparse becomes CSR, dense stays dense).
    """
    order = np.argsort(row_key)
    unsort = np.empty_like(order)
    unsort[order] = np.arange(order.size)
    X_rows = mat[row_key[order], :]
    if sp.issparse(X_rows):
        return X_rows.tocsr()[unsort, :]
    return X_rows[unsort, :]


class AnnDataCollection:
    """Wrapper representing a merged collection of AnnData files.
    Example:
        ds = AnnDataCollection.from_files(["a.h5ad", "b.h5ad"], out_path="merged.h5ad")
        view = ds[[True]*100, ["GeneA", "GeneB"]]
        adata = view.to_memory()
    """

    def __init__(
        self,
        adata_paths: List[str],
        adata: Optional[anndata.AnnData] = None,
        source_info: Optional[List[Dict]] = None,
    ):
        """Initialize an `AnnDataCollection` wrapper.
        Parameters
        - `adata_paths`: list of filesystem paths to the individual `.h5ad` files.
        - `adata`: optional `anndata.AnnData` that contains merged metadata
            (typically `obs`, `var`, and `obsm`) for the combined dataset. When
            present this `AnnData` usually contains an empty or placeholder
            `X` matrix because the heavy expression matrices remain stored in
            the original files and are loaded on demand by `AnnSubset.to_memory()`.
            If `ann` is `None`, the wrapper does not have merged metadata yet.
        - `source_info`: optional list of dictionaries with per-source
            metadata (e.g. path, n_obs, var_names). Populated by
            `AnnDataCollection.from_files()` when creating a merged dataset.
        """
        self.adata_paths = list(adata_paths)
        self.adata = adata
        self.source_info = source_info or []

    @classmethod
    def from_files(
        cls,
        paths: Sequence[str],
        out_path: Optional[str] = None,
        metadata_path: Optional[str] = None,
        join: str = "inner",
    ) -> "AnnDataCollection":
        """Create a merged AnnDataCollection from existing `.h5ad` files.
        This will merge `obs` (stacked) and `var` across files according to
        `join`. `X` is not merged; the saved on-disk AnnData contains an empty
        sparse matrix of shape (n_obs_total, n_vars_total). The merged var order
        is what a subsequent `write_h5ad` (or `to_memory`) reindexes every
        source's columns to.

        Parameters
        ----------
        join : str, optional
            How to merge var_names across sources. ``"inner"`` (default) keeps
            the intersection (columns present in every file); ``"outer"`` keeps
            the union (columns present in any file, missing ones filled with 0
            when the matrix is later streamed). Both preserve the first file's
            var order (outer appends any new vars in first-seen order). Sources
            with identical vars behave the same under either.
        """
        # `paths` should be an iterable of paths to `.h5ad` files, all .h5ad files must
        # have the same vars. When
        # `out_path` is provided the merged AnnData (metadata-only) will be
        # written to that file. The returned `AnnDataCollection` keeps the list of
        # original file paths and stores per-source info in `source_info`.
        # Allow `paths` to be either:
        # - a single wildcard string: "/path/*.h5ad"
        # - a sequence of explicit paths
        # - a sequence with one wildcard string
        if isinstance(paths, str):
            paths = sorted(glob.glob(os.path.expanduser(paths)))
        else:
            # convert to list of str
            paths = [str(p) for p in paths]

        if not paths:
            raise FileNotFoundError("No .h5ad files found for given paths")
        if join not in ("inner", "outer"):
            raise ValueError(f"join must be 'inner' or 'outer', got {join!r}")
        obs_list = []
        # preserve var order as seen across files for faster subsetting later
        # var_list: List[str] = []
        # _seen_vars = set()
        per_source_info: List[Dict[str, Any]] = []
        var_names = None

        def _read_source_meta(i_p):
            """Read obs metadata from a single source file."""
            i, p = i_p
            a = anndata.read_h5ad(p, backed="r")
            try:
                obs = a.obs.copy()
                obs["_orig_obs_name"] = a.obs_names.astype(str)
                obs["_orig_obs_pos"] = np.arange(a.n_obs, dtype=np.int64)
                obs["_source_idx"] = i
                info = {
                    "path": os.path.abspath(p),
                    "n_obs": a.n_obs,
                    "n_vars": a.n_vars,
                }
                source_var_names = list(a.var_names)
                return i, obs, info, source_var_names
            finally:
                _close_backed(a)

        # First pass: gather obs frames and per-source info.
        # Read in backed mode so we don't load X into memory.
        # Parallelize across source files for faster metadata gathering.
        n_workers = min(8, max(1, len(paths)))
        results_by_idx: Dict[int, Any] = {}
        with ThreadPoolExecutor(max_workers=n_workers) as exe:
            for i, obs, info, src_vars in exe.map(_read_source_meta, enumerate(paths)):
                results_by_idx[i] = (obs, info, src_vars)
        # Process results in order to preserve deterministic obs ordering
        source_var_lists: List[List[str]] = []
        for i in range(len(paths)):
            obs, info, src_vars = results_by_idx[i]
            obs_list.append(obs)
            per_source_info.append(info)
            source_var_lists.append(list(src_vars))

        # Merge var_names across sources according to `join`. Both modes
        # preserve the first file's var order; the streaming write/read reindex
        # each source's columns to this order (filling absent ones with 0 for
        # 'outer'). Identical-var sources give the same result either way.
        if join == "inner":
            common = set(source_var_lists[0])
            for sv in source_var_lists[1:]:
                common &= set(sv)
            var_names = [v for v in source_var_lists[0] if v in common]
            if not var_names:
                raise ValueError(
                    "join='inner' produced an empty var intersection across the "
                    "input files; use join='outer' or check the var_names."
                )
        else:  # outer: union, first-seen order across files
            seen = set()
            var_names = []
            for sv in source_var_lists:
                for v in sv:
                    if v not in seen:
                        seen.add(v)
                        var_names.append(v)

        # Concatenate obs
        merged_obs = pd.concat(obs_list)
        merged_obs.index = merged_obs.index.map(str)

        # If metadata_path provided, read metadata and filter merged_obs to only those cells.
        if metadata_path is not None:
            sep = (
                "\t"
                if metadata_path.endswith(".tsv")
                or metadata_path.endswith(".tsv.gz")
                or metadata_path.endswith(".txt")
                else ","
            )
            metadata_path = os.path.expanduser(metadata_path)
            # low_memory=False reads the whole file at once so column dtypes are
            # inferred consistently (avoids the chunked-inference DtypeWarning on
            # columns that look numeric in some rows and string in others).
            meta = pd.read_csv(metadata_path, index_col=0, sep=sep, low_memory=False)
            meta_ids = set(meta.index.astype(str).tolist())
            # keep_ids= list(set(merged_obs.index.astype(str)).intersection(meta_ids))
            keep_ids = [
                cell
                for cell in merged_obs.index.astype(str).tolist()
                if cell in meta_ids
            ]
            merged_obs = merged_obs.loc[keep_ids]
            for col in meta.columns:
                if col not in merged_obs.columns:
                    merged_obs[col] = meta.loc[merged_obs.index.tolist(), col].tolist()

        # Build merged var DataFrame (union of var names)
        merged_var = pd.DataFrame(index=pd.Index(var_names))  # var_list
        # recompute total observations after optional metadata filtering
        # Create an empty sparse X with proper shape (we don't merge expression matrices here)
        X_empty = sp.csr_matrix(
            (merged_obs.shape[0], len(merged_var)), dtype=np.float32
        )
        adata = anndata.AnnData(X_empty, obs=merged_obs, var=merged_var)
        # Note: we intentionally do not merge `obsm` arrays. If per-source
        # `obsm` data is needed, read individual files directly from
        # `individual_adata_paths` stored in `merged.uns`.

        # store metadata (short keys)
        adata.uns["src_paths"] = [os.path.abspath(p) for p in paths]
        # JSON-serialize per-source info so HDF5 can store it safely
        try:
            adata.uns["src_info"] = json.dumps(per_source_info)
        except Exception:
            # fallback: store as a list of strings
            adata.uns["src_info"] = [str(x) for x in per_source_info]
        if metadata_path is not None:
            adata.uns["metadata_path"] = os.path.abspath(metadata_path)
        adata.uns["join"] = join
        instance = cls(list(paths), adata=adata, source_info=per_source_info)
        # optionally save to out_path
        if out_path:
            adata.write_h5ad(out_path)
        return instance

    @classmethod
    def read(cls, path: str) -> "AnnDataCollection":
        """Load a merged `.h5ad` file created by `from_files` and return an
        `AnnDataCollection`.

        Expects `src_paths` (or `individual_adata_paths`) in `adata.uns` and
        `src_info` which may be JSON-serialized or a list.
        """
        path = os.path.expanduser(path)
        adata = anndata.read_h5ad(path)

        # Basic validation: ensure this AnnData was produced by `from_files`.
        has_src_uns = ("src_paths" in adata.uns) or (
            "individual_adata_paths" in adata.uns
        )
        has_obs_markers = ("_source_idx" in adata.obs.columns) or (
            "_orig_obs_name" in adata.obs.columns
        )
        if not (has_src_uns and has_obs_markers):
            raise ValueError(
                f"File {path} does not appear to be an AnnDataCollection (missing src_paths/src_info or obs markers)"
            )

        # Avoid relying on truthiness of containers (e.g., numpy arrays)
        sp_raw = adata.uns.get("src_paths", None)
        if sp_raw is None:
            sp_raw = adata.uns.get("individual_adata_paths", None)
        if sp_raw is None:
            src_paths = []
        else:
            # normalize to a list of strings regardless of input type
            if isinstance(sp_raw, (list, tuple, pd.Index)):
                src_paths = [str(p) for p in sp_raw]
            elif isinstance(sp_raw, np.ndarray):
                src_paths = [str(p) for p in sp_raw.tolist()]
            else:
                src_paths = [str(sp_raw)]

        raw_info = adata.uns.get("src_info", None)
        source_info: List[Dict[str, Any]] = []
        if raw_info is None:
            source_info = []
        elif isinstance(raw_info, str):
            try:
                source_info = json.loads(raw_info)
            except Exception:
                source_info = [raw_info]
        elif isinstance(raw_info, list):
            source_info = raw_info
        else:
            source_info = [raw_info]

        return cls(src_paths, adata=adata, source_info=source_info)

    def write_h5ad(
        self,
        path: str,
        thread: int = 8,
        layer: Optional[str] = None,
        dense: Optional[bool] = None,
        compression: Optional[str] = "gzip",
        compression_opts: Optional[int] = None,
        block_size: Optional[int] = None,
        max_block_gb: float = 2.0,
    ) -> None:
        """Write this collection's real `X` to a consolidated `.h5ad`.

        Memory-efficient: writes the obs/var/uns/obsm metadata skeleton first
        (with an empty `X`), then rebuilds the `X` group in place by reading a
        bounded block of rows at a time from `adata_paths` and appending to
        (dense in-place or resizable CSR) HDF5 datasets. Each source file is
        streamed in row blocks (NOT loaded whole), so peak memory is ~`thread`
        blocks and stays bounded no matter how large a single source file is.
        The result reads back as a plain `AnnData`.

        Parameters
        ----------
        path : str
            Destination `.h5ad` path.
        thread : int, optional
            Number of row blocks read concurrently (bounded prefetch). Peak
            memory is ~``thread * max_block_gb``. Default 8.
        layer : str or None, optional
            Name of the source layer to read from each `.h5ad` file. When
            `None` (default) the main `X` matrix is used; otherwise the named
            entry in `layers` is read. The written matrix is always stored as
            the output file's `X`.
        dense : bool or None, optional
            Storage format for the written `X`. When ``True`` the matrix is
            written as a dense HDF5 array; when ``False`` as a streamed CSR
            group. When ``None`` (default) the format is auto-detected from the
            first source file's `X` (or `layer`) encoding, so a dense source
            stays dense and a sparse source stays sparse.
        compression : str or None, optional
            HDF5 compression filter applied to both the metadata skeleton and
            the streamed `X` datasets. One of ``"gzip"`` (default), ``"lzf"``,
            or ``None`` to disable compression. ``"gzip"`` gives the smallest
            files (portable), ``"lzf"`` is faster but compresses less.
        compression_opts : int or None, optional
            Filter tuning passed to h5py. For ``"gzip"`` this is the level 0-9
            (defaults to 4 when ``None``); ignored for ``"lzf"``.
        block_size : int or None, optional
            Number of rows (cells) read per block from each source. When
            ``None`` (default) it is derived from ``max_block_gb`` and the
            matrix width so one dense block stays within the budget; the wider
            the matrix, the fewer rows per block. Set an explicit value to
            override.
        max_block_gb : float, optional
            Approximate memory budget (in GB) for a single row block when
            ``block_size`` is None. Peak memory is ~``thread * max_block_gb``.
            Default 2.0. Lower it (or ``thread``) if memory is tight; raise it
            for fewer, larger IO passes on narrow matrices.
        """
        import h5py

        if self.adata is None:
            raise ValueError(
                "AnnDataCollection has no merged metadata (`adata` is None); "
                "nothing to write."
            )
        path = os.path.expanduser(path)

        merged = self.adata
        obs = merged.obs
        n_obs = obs.shape[0]
        out_var_names = list(map(str, merged.var_names))
        n_vars = len(out_var_names)
        paths = self.adata_paths

        if "_source_idx" not in obs.columns:
            raise ValueError(
                "AnnDataCollection obs is missing '_source_idx'; cannot locate "
                "source rows for streaming write."
            )
        src_indices = obs["_source_idx"].values.astype(int)
        has_pos = "_orig_obs_pos" in obs.columns
        if has_pos:
            orig_pos = obs["_orig_obs_pos"].values.astype(int)
        elif "_orig_obs_name" in obs.columns:
            orig_name = obs["_orig_obs_name"].to_numpy()
        else:
            raise ValueError(
                "AnnDataCollection obs is missing both '_orig_obs_pos' and "
                "'_orig_obs_name'; cannot locate source rows."
            )

        # Rows per block: bound one dense block to ~max_block_gb so a very wide
        # matrix (e.g. genome-wide 5kb, hundreds of thousands of cols) uses few
        # rows per read while a narrow one uses many. float32 output -> 4 bytes.
        if block_size is not None:
            rows_per_block = max(1, int(block_size))
        else:
            budget_bytes = max(1, int(float(max_block_gb) * 1e9))
            rows_per_block = max(1, budget_bytes // max(1, n_vars * 4))

        # Build contiguous same-source runs over the final obs order, split into
        # bounded row blocks so each read (and the write_block that follows) only
        # materializes `rows_per_block` rows, not a whole source file. Kept in
        # increasing-start order so the CSR append path stays row-ordered.
        runs = []  # (start, end, source_idx, row_key)
        i = 0
        while i < n_obs:
            s = src_indices[i]
            j = i + 1
            while j < n_obs and src_indices[j] == s:
                j += 1
            for bstart in range(i, j, rows_per_block):
                bend = min(bstart + rows_per_block, j)
                row_key = (orig_pos[bstart:bend] if has_pos
                           else orig_name[bstart:bend])
                runs.append((bstart, bend, int(s), row_key))
            i = j

        # Precompute each source's column map ONCE (cheap h5py open, ~1 ms)
        # so the many per-block reads don't re-parse the source's (possibly
        # huge, e.g. 460k) var index every time via anndata.read_h5ad (~0.5 s).
        col_maps = {}
        for s in sorted({r[2] for r in runs}):
            try:
                col_maps[s] = _source_col_map(paths[s], out_var_names, layer)
            except Exception:
                col_maps[s] = None  # reader falls back to the anndata path

        # Normalize compression settings shared by the metadata skeleton and
        # the streamed X datasets. gzip takes an int level via compression_opts;
        # lzf takes none. Build kwargs once for the h5py create_dataset calls.
        if compression == "gzip" and compression_opts is None:
            compression_opts = 4
        ds_comp_kwargs = {}
        if compression is not None:
            ds_comp_kwargs["compression"] = compression
            if compression == "gzip":
                ds_comp_kwargs["compression_opts"] = compression_opts

        # 1) Write metadata skeleton (empty X) via anndata so obs/var/uns/obsm
        #    are encoded correctly. Drop collection-identifying uns keys so the
        #    result reads back as a plain AnnData.
        meta = merged.copy()
        meta.X = sp.csr_matrix((n_obs, n_vars), dtype=np.float32)
        for k in ("src_paths", "individual_adata_paths", "src_info", "metadata_path"):
            meta.uns.pop(k, None)
        if compression is not None:
            meta.write_h5ad(path, compression=compression,
                            compression_opts=compression_opts)
        else:
            meta.write_h5ad(path)

        # 2) Replace the X group/dataset with a streamed matrix. Output density
        #    follows `dense` when given, else the first source's encoding.
        if dense is None:
            dense = (
                _detect_x_encoding(paths[runs[0][2]], layer) == "dense"
                if runs
                else False
            )

        with h5py.File(path, "a") as f:
            if "X" in f:
                del f["X"]
            if dense:
                # Dense output: shape is known up front, so fill row blocks in
                # place (no resize needed). Peak memory stays bounded to a few
                # source blocks.
                X_ds = f.create_dataset(
                    "X", shape=(n_obs, n_vars), dtype=np.float32, chunks=True,
                    **ds_comp_kwargs
                )
                X_ds.attrs["encoding-type"] = "array"
                X_ds.attrs["encoding-version"] = "0.2.0"

                def _write_block(start, end, block):
                    X_ds[start:end, :] = block

                _stream_source_blocks(
                    runs, paths, out_var_names, layer, thread, True, _write_block,
                    col_maps=col_maps
                )
            else:
                # Sparse output: append each block's data/indices to resizable
                # datasets and accumulate the shared indptr.
                indptr = np.zeros(n_obs + 1, dtype=np.int64)
                g = f.create_group("X")
                g.attrs["encoding-type"] = "csr_matrix"
                g.attrs["encoding-version"] = "0.1.0"
                g.attrs["shape"] = np.array([n_obs, n_vars], dtype=np.int64)
                data_ds = g.create_dataset(
                    "data", shape=(0,), maxshape=(None,), dtype=np.float32,
                    chunks=True, **ds_comp_kwargs
                )
                indices_ds = g.create_dataset(
                    "indices", shape=(0,), maxshape=(None,), dtype=np.int32,
                    chunks=True, **ds_comp_kwargs
                )

                nnz = 0

                def _write_block(start, end, block):
                    nonlocal nnz
                    bd = block.data
                    bi = block.indices.astype(np.int32)
                    m = bd.shape[0]
                    data_ds.resize((nnz + m,))
                    data_ds[nnz : nnz + m] = bd
                    indices_ds.resize((nnz + m,))
                    indices_ds[nnz : nnz + m] = bi
                    indptr[start + 1 : end + 1] = block.indptr[1:] + nnz
                    nnz += m

                _stream_source_blocks(
                    runs, paths, out_var_names, layer, thread, False, _write_block,
                    col_maps=col_maps
                )

                g.create_dataset("indptr", data=indptr, dtype=np.int64,
                                 **ds_comp_kwargs)

    def __len__(self) -> int:
        """Return the number of observations (cells) in the collection."""
        return self.adata.n_obs if self.adata is not None else 0

    def __getitem__(self, key: Tuple[Any, Any]) -> "AnnDataView":
        """Support `adatas[cells, genes]` slicing.

        `cells` and `genes` may be:
          - slice, list of booleans, list of indices, list of labels
        Returns an `AnnDataView` with a `to_memory()` method that will load X.
        """
        if not isinstance(key, tuple) or len(key) != 2:
            raise IndexError("Indexing must be adatas[cells, genes]")

        cells, genes = key
        obs_idx = self._parse_obs_indexer(cells)
        var_names = self._parse_var_indexer(genes)
        return AnnDataView(self, obs_idx, var_names)

    def update_paths(
        self,
        mapping: Optional[Dict[str, str]] = None,
        search_dirs: Optional[Sequence[str]] = None,
        recursive: bool = True,
    ) -> Dict[str, Optional[str]]:
        """Update stored source paths after files have been moved.

        Parameters
        ----------
        mapping : dict or None, optional
            Dict mapping old_path -> new_path.  If provided, entries in
            ``adata_paths`` will be replaced according to this map.
        search_dirs : sequence of str or None, optional
            Directories to search for moved files by matching basenames
            (first match wins).  Used when *mapping* is not provided or
            doesn't contain an entry for a given source.
        recursive : bool, optional
            If ``True``, search directories recursively when using
            *search_dirs*.  Default is ``True``.

        Returns
        -------
        dict
            Mapping of original absolute path to new absolute path
            (or ``None`` if not found).
        """
        result: Dict[str, Optional[str]] = {}

        norm_map: Dict[str, str] = {}
        if mapping:
            for k, v in mapping.items():
                try:
                    norm_map[os.path.abspath(k)] = os.path.abspath(v)
                except Exception:
                    norm_map[str(k)] = str(v)

        search_dirs = [os.path.expanduser(d) for d in (search_dirs or [])]
        for i, orig in enumerate(list(self.adata_paths)):
            orig_abs = os.path.abspath(orig)
            new_path: Optional[str] = None

            # 1) explicit mapping provided by user
            if orig_abs in norm_map:
                cand = norm_map[orig_abs]
                if os.path.exists(cand):
                    new_path = cand

            # 2) search by basename in provided search_dirs
            if new_path is None and search_dirs:
                basename = os.path.basename(orig_abs)
                for d in search_dirs:
                    if recursive:
                        pattern = os.path.join(d, "**", basename)
                    else:
                        pattern = os.path.join(d, basename)
                    matches = sorted(glob.glob(pattern, recursive=recursive))
                    if matches:
                        new_path = os.path.abspath(matches[0])
                        break

            # 3) if original still exists, keep it
            if new_path is None and os.path.exists(orig_abs):
                new_path = orig_abs

            result[orig_abs] = new_path

            if new_path:
                # update internal list
                self.adata_paths[i] = new_path
                # update adata.uns if present (accept old key as fallback)
                if self.adata is not None:
                    try:
                        paths_list = list(
                            self.adata.uns.get(
                                "src_paths",
                                self.adata.uns.get("individual_adata_paths", []),
                            )
                        )
                        if i < len(paths_list):
                            paths_list[i] = new_path
                        else:
                            # ensure list is long enough
                            while len(paths_list) < i:
                                paths_list.append("")
                            paths_list.append(new_path)
                        self.adata.uns["src_paths"] = [
                            os.path.abspath(p) for p in paths_list
                        ]
                    except Exception:
                        # best-effort; don't raise on metadata update
                        pass
                # update source_info if present
                if i < len(self.source_info):
                    try:
                        self.source_info[i]["path"] = new_path
                    except Exception:
                        pass

        # re-serialize per-source info into uns for safe HDF5 storage
        if self.adata is not None and self.source_info:
            try:
                self.adata.uns["src_info"] = json.dumps(self.source_info)
            except Exception:
                self.adata.uns["src_info"] = [str(x) for x in self.source_info]

        return result

    def _parse_obs_indexer(self, idx) -> np.ndarray:
        """
        Convert an observation indexer to integer position indices.

        Parameters
        ----------
        idx : None, slice, list, np.ndarray, or pd.Series
            - None: select all observations.
            - slice: standard Python slice.
            - bool array: select where True.
            - int array: select by position.
            - str array: select by obs_names (cell IDs).

        Returns
        -------
        np.ndarray
            Integer position indices into the observation axis.
        """
        n = self.adata.n_obs
        if idx is None:
            return np.arange(n)
        if isinstance(idx, slice):
            return np.arange(n)[idx]
        if isinstance(idx, (list, np.ndarray, pd.Series)):
            arr = np.asarray(idx)
            if arr.dtype == bool:
                return np.nonzero(arr)[0]  # get indexes where True
            # If integer dtype, treat as indices. Otherwise treat as labels.
            if np.issubdtype(arr.dtype, np.integer):
                return arr.astype(int)
            # treat as obs names (strings or object)
            return np.array(
                self.adata.obs_names.get_indexer(arr.astype(str)), dtype=int
            )
        raise IndexError("Unsupported obs indexer")

    def _parse_var_indexer(self, idx) -> List[str]:
        """
        Convert a variable indexer to a list of variable names.

        Parameters
        ----------
        idx : None, slice, list, np.ndarray, pd.Series, or str
            - None: select all variables.
            - slice: standard Python slice on var_names.
            - bool array: select where True.
            - int array: select by position.
            - str array/list: select by variable names directly.
            - str: a single variable name.

        Returns
        -------
        list of str
            List of selected variable names.
        """
        if idx is None:  # idx==slice(None)
            return list(self.adata.var_names)
        if isinstance(idx, slice):
            return list(self.adata.var_names[idx])
        if isinstance(idx, (list, np.ndarray, pd.Series)):
            arr = np.asarray(idx)
            if arr.dtype == bool:
                return list(np.asarray(self.adata.var_names)[arr])
            # If integer dtype, interpret as indices into var_names
            if np.issubdtype(arr.dtype, np.integer):
                return list(np.asarray(self.adata.var_names)[arr.astype(int)])
            # Otherwise treat as labels (string/object/unicode)
            return list(arr.astype(str))
        if isinstance(idx, str):
            return [idx]
        raise IndexError("Unsupported var indexer")


def is_anndatacollection(path: str) -> bool:
    """Return True if `path` points to an AnnData file produced by
    `AnnDataCollection.from_files`.
    Heuristic: the AnnData must have `uns` entries `src_paths` or
    `individual_adata_paths` and its `obs` must contain either
    `_source_idx` or `_orig_obs_name` markers.
    """
    if path is None:
        return False
    path = os.path.expanduser(path)
    try:
        ad = anndata.read_h5ad(path, backed="r")
    except Exception:
        return False
    try:
        has_src_uns = ("src_paths" in ad.uns) or ("individual_adata_paths" in ad.uns)
        # obs may be empty; check columns safely
        has_obs_markers = False
        try:
            has_obs_markers = ("_source_idx" in ad.obs.columns) or (
                "_orig_obs_name" in ad.obs.columns
            )
        except Exception:
            has_obs_markers = False
        return bool(has_src_uns and has_obs_markers)
    finally:
        _close_backed(ad)


class AnnDataView:
    """Helper representing a (cells, genes) subset of an `AnnDataCollection`.
    Call `to_memory()` to read and assemble the actual `AnnData` with `X`.
    """

    def __init__(
        self,
        dataset: AnnDataCollection,
        obs_idx: np.ndarray = None,
        var_names: List[str] = None,
    ):
        """
        Initialize an AnnDataView.

        Parameters
        ----------
        dataset : AnnDataCollection
            Parent collection to read data from.
        obs_idx : np.ndarray, optional
            Integer indices of selected observations. If None, all
            observations are selected.
        var_names : list of str, optional
            Names of selected variables. If None, all variables are
            selected.
        """
        self.dataset = dataset
        self.obs_idx = (
            np.asarray(obs_idx, dtype=int) if obs_idx is not None else None
        )  # np.arange(dataset.adata.n_obs)
        self.var_names = (
            list(var_names) if var_names is not None else None
        )  # list(dataset.adata.var_names)

    def to_memory(self, thread=8, layer=None) -> anndata.AnnData:
        """Load selected `X` chunks from underlying files and assemble AnnData.
        Returns a new `anndata.AnnData` with concatenated `X` for the selected
        cells and genes.

        Parameters
        ----------
        thread : int, optional
            Number of source files read concurrently. Default 8.
        layer : str or None, optional
            Name of the source layer to read from each `.h5ad` file. When
            `None` (default) the main `X` matrix is used; otherwise the named
            entry in `layers` is read. The extracted matrix is always placed in
            the returned AnnData's `X`.
        """
        ds = self.dataset
        merged = ds.adata

        # Build target obs and var
        if self.obs_idx is not None:
            sel_obs = merged.obs.iloc[self.obs_idx].copy()
        else:
            sel_obs = merged.obs.copy()
        if self.var_names is not None:
            sel_var = merged.var.loc[self.var_names].copy()
        else:
            sel_var = merged.var.copy()

        n_rows = len(sel_obs)
        n_cols = len(sel_var)

        # Map merged obs rows to sources and cache locals for speed
        src_indices = sel_obs["_source_idx"].values.astype(int)
        paths = ds.adata_paths

        # `_orig_obs_pos` and `_orig_obs_name` are per-cell markers written by
        # `from_files`: they record where each merged row came from in its
        # source file. `_orig_obs_pos` is the integer row index (0-based
        # position) in the source's obs axis, which allows fast positional HDF5
        # reads. `_orig_obs_name` is the original obs_name (cell ID) string,
        # used as a slower label-based fallback when positions are unavailable.
        # Check for integer position column (faster HDF5 access)
        has_pos = "_orig_obs_pos" in sel_obs.columns
        has_name = "_orig_obs_name" in sel_obs.columns
        if has_pos:
            sel_orig_pos = sel_obs["_orig_obs_pos"].values.astype(int)
        if has_name:
            sel_orig_obs = sel_obs["_orig_obs_name"].to_numpy()
        if not (has_pos or has_name):
            raise ValueError(
                "sel_obs is missing both '_orig_obs_pos' and '_orig_obs_name'; "
                "cannot locate source rows."
            )

        # Build column mapping: positions in source var axis → positions in output var axis.
        # All sources share the same var_names, so we compute this once.
        var_names_arr = sel_var.index.astype(str).tolist()
        merged_var_map = {v: idx for idx, v in enumerate(merged.var_names)}  # type: ignore
        src_gene_pos_all = np.array(
            [merged_var_map.get(v, -1) for v in var_names_arr], dtype=int
        )
        present_mask = src_gene_pos_all >= 0
        if not np.any(present_mask):
            return None
        # Positions in the *output* matrix for each extracted gene
        src_col_pos_in_output = np.nonzero(present_mask)[0] # such as array([0, 2, 3])
        # Positions in each *source* file's var axis
        src_gene_pos = src_gene_pos_all[present_mask]

        # Fast path: when every source column is selected in natural order the
        # CSR column-gather below is a pure (and costly) copy — skip it.
        n_source_vars = len(merged.var_names)
        full_cols = src_gene_pos.size == n_source_vars and np.array_equal(
            src_gene_pos, np.arange(n_source_vars)
        )

        # Pre-compute per-source work items (avoid submitting idle tasks)
        unique_sources = np.unique(src_indices)
        work_items = []  # (src_idx, path, global_rows, row_positions_or_names)
        for src_idx in unique_sources:
            if src_idx < 0 or src_idx >= len(paths):
                continue
            mask = src_indices == src_idx
            global_rows = np.nonzero(mask)[0]
            if has_pos:
                row_key = sel_orig_pos[mask]  # integer positions — fast
            else:
                row_key = sel_orig_obs[mask]  # obs names — fallback
            work_items.append((int(src_idx), paths[src_idx], global_rows, row_key))

        # Validate the requested layer once up front. Per-source read errors are
        # swallowed below (partial-result tolerance), so without this a mistyped
        # layer name would silently yield an empty matrix.
        if layer is not None and work_items:
            _vp = work_items[0][1]
            _va = anndata.read_h5ad(_vp, backed="r")
            try:
                if layer not in _va.layers:
                    raise ValueError(
                        f"Layer '{layer}' not found in source file {_vp}. "
                        f"Available layers: {list(_va.layers.keys())}"
                    )
            finally:
                _close_backed(_va)

        def _process_source(item):
            """Read and extract X data from a single source h5ad file."""
            src_idx, path, global_rows, row_key = item
            a = anndata.read_h5ad(path, backed="r")
            try:
                # Read selected rows from disk
                if np.issubdtype(row_key.dtype, np.integer):
                    X_rows = _read_rows_in_order(_matrix_row_source(a, layer), row_key)
                else:
                    sub = a[row_key, :]
                    X_rows = sub.X if layer is None else sub.layers[layer]

                # Ensure CSR for fast column slicing
                if not sp.issparse(X_rows):
                    X_rows = sp.csr_matrix(X_rows)
                elif not sp.isspmatrix_csr(X_rows):
                    X_rows = X_rows.tocsr()

                # Extract only the needed columns (skip the copy for full slices)
                X_sub = X_rows if full_cols else X_rows[:, src_gene_pos]
                X_coo = sp.coo_matrix(X_sub)

                if X_coo.nnz == 0:
                    return None

                rows_mapped = global_rows[X_coo.row].astype(np.int32)
                cols_mapped = src_col_pos_in_output[X_coo.col].astype(np.int32)
                data_arr = X_coo.data.astype(np.float32)

                return rows_mapped, cols_mapped, data_arr
            finally:
                _close_backed(a)

        # Process only sources that have selected cells, in parallel
        rows_all: List[np.ndarray] = []
        cols_all: List[np.ndarray] = []
        data_all: List[np.ndarray] = []

        max_workers = min(thread, max(1, len(work_items)))
        if len(work_items) == 0:
            pass  # no work to do
        elif len(work_items) == 1:
            # Skip ThreadPoolExecutor overhead for single source
            res = _process_source(work_items[0])
            if res is not None:
                rows_all.append(res[0])
                cols_all.append(res[1])
                data_all.append(res[2])
        else:
            with ThreadPoolExecutor(max_workers=max_workers) as exe:
                futures = {
                    exe.submit(_process_source, item): item for item in work_items
                }
                for fut in as_completed(futures):
                    try:
                        res = fut.result()
                    except Exception:
                        logger.error("Error processing source %s", futures[fut][1])
                        continue
                    if res is not None:
                        rows_all.append(res[0])
                        cols_all.append(res[1])
                        data_all.append(res[2])

        if len(data_all) == 0:
            X_final = sp.csr_matrix((n_rows, n_cols), dtype=np.float32)
        else:
            rows_cat = np.concatenate(rows_all)
            cols_cat = np.concatenate(cols_all)
            data_cat = np.concatenate(data_all)
            X_final = sp.coo_matrix(
                (data_cat, (rows_cat, cols_cat)), shape=(n_rows, n_cols)
            ).tocsr()

        out = anndata.AnnData(X_final, obs=sel_obs, var=sel_var)
        return out



def read_h5ad(path: str, **kwargs):
    """Read a `.h5ad` file, dispatching on its type.

    If `path` points to a merged file produced by
    `AnnDataCollection.from_files` (detected via `is_anndatacollection`), it is
    loaded as an `AnnDataCollection`. Otherwise it is read as a plain
    `anndata.AnnData` via `anndata.read_h5ad`.

    Parameters
    ----------
    path : str
        Path to the `.h5ad` file.
    **kwargs
        Extra keyword arguments forwarded to `anndata.read_h5ad` when the file
        is a plain `AnnData` (e.g. `backed="r"`). Ignored for
        `AnnDataCollection` files.

    Returns
    -------
    AnnDataCollection or anndata.AnnData
    """
    path = os.path.expanduser(path)
    if is_anndatacollection(path):
        return AnnDataCollection.read(path)
    return anndata.read_h5ad(path, **kwargs)


def _detect_x_encoding(path: str, layer: Optional[str] = None) -> str:
    """Return ``"dense"`` or ``"sparse"`` for a source file's X (or named layer).

    Sparse matrices are stored as an HDF5 group (csr/csc); dense arrays as an
    HDF5 dataset. Falls back to ``"sparse"`` if the element cannot be inspected.
    """
    a = anndata.read_h5ad(path, backed="r")
    try:
        elem = a.file["X"] if layer is None else a.file["layers"][layer]
        return "dense" if isinstance(elem, h5py.Dataset) else "sparse"
    except Exception:
        return "sparse"
    finally:
        _close_backed(a)


def _read_h5_index_names(grp) -> List[str]:
    """Read the index (row/col names) of an anndata obs/var h5py group cheaply."""
    idx = grp.attrs.get("_index", "_index")
    if isinstance(idx, bytes):
        idx = idx.decode()
    vals = grp[idx][:]
    return [v.decode() if isinstance(v, (bytes, bytearray)) else str(v) for v in vals]


def _source_col_map(path: str, out_var_names, layer):
    """Precompute a source file's column mapping to the output var order.

    Opens the file once via ``h5py`` (cheap: ~1 ms, and crucially does NOT parse
    the whole obs/var DataFrames the way ``anndata.read_h5ad`` does — that costs
    ~0.5 s per open on a 460k-var file). Returns ``(col_pos, present, identity,
    is_dense)`` where ``col_pos[k]`` is the source column index for output var
    ``k`` (or -1 if absent), ``present`` is the bool mask, ``identity`` means the
    source vars already equal the output order (no gather needed), and
    ``is_dense`` flags a dense on-disk matrix.
    """
    with h5py.File(path, "r") as f:
        src_var = _read_h5_index_names(f["var"])
        elem = f["X"] if layer is None else f["layers"][layer]
        is_dense = isinstance(elem, h5py.Dataset)
    pos_map = {v: i for i, v in enumerate(src_var)}
    col_pos = np.array([pos_map.get(v, -1) for v in out_var_names], dtype=np.int64)
    present = col_pos >= 0
    identity = src_var == list(out_var_names)
    return col_pos, present, identity, is_dense


def _read_block_fast(path, row_pos, col_map, layer=None, as_dense=False):
    """Read integer-position ``row_pos`` rows of a source via a cheap h5py open.

    Uses the precomputed ``col_map`` (see :func:`_source_col_map`) so no per-block
    var-name reparse is needed; opening with h5py is ~1 ms, so streaming a file in
    many small blocks stays fast. Only one block (``len(row_pos) x n_out``) is
    held in memory. Returns a dense float32 ndarray when ``as_dense`` else CSR.
    """
    col_pos, present, identity, is_dense = col_map
    row_pos = np.asarray(row_pos)
    order = np.argsort(row_pos, kind="stable")
    unsort = np.empty_like(order)
    unsort[order] = np.arange(order.size)
    sorted_rows = row_pos[order]
    with h5py.File(path, "r") as f:
        elem = f["X"] if layer is None else f["layers"][layer]
        if is_dense:
            X = np.asarray(elem[sorted_rows, :])[unsort, :]
        else:
            if _sparse_dataset is None:
                raise RuntimeError("sparse_dataset unavailable; use fallback reader")
            X = _sparse_dataset(elem)[sorted_rows, :]
            X = (X.tocsr() if sp.issparse(X) else sp.csr_matrix(X))[unsort, :]
    n_out = int(col_pos.shape[0])
    if identity:
        out = X
    elif present.all():
        out = X[:, col_pos]
    else:
        cols_in = np.nonzero(present)[0]
        if sp.issparse(X):
            out = sp.lil_matrix((X.shape[0], n_out), dtype=np.float32)
            out[:, cols_in] = X[:, col_pos[present]]
            out = out.tocsr()
        else:
            out = np.zeros((X.shape[0], n_out), dtype=X.dtype)
            out[:, cols_in] = X[:, col_pos[present]]
    if as_dense:
        return np.asarray(out.toarray() if sp.issparse(out) else out, dtype=np.float32)
    if not sp.issparse(out):
        out = sp.csr_matrix(out)
    return out.astype(np.float32, copy=False)


def _stream_source_blocks(
    runs,
    paths,
    out_var_names,
    layer,
    thread,
    as_dense,
    write_block,
    col_maps=None,
) -> None:
    """Read source blocks in run order and hand each to ``write_block``.

    Uses bounded prefetch (at most ``thread`` blocks in flight) so blocks are
    consumed sequentially while reads overlap. Each block passed to
    ``write_block(start, end, block)`` is a CSR matrix when ``as_dense`` is
    False, otherwise a dense float32 ndarray. When ``col_maps`` (a
    ``{source_idx: col_map}`` dict from :func:`_source_col_map`) is given and the
    row key is integer positions, the fast h5py reader is used (no per-block
    var-name reparse); otherwise it falls back to the robust anndata reader.
    """
    from collections import deque

    def _read(s, rk):
        rk = np.asarray(rk)
        if (
            col_maps is not None
            and col_maps.get(s) is not None
            and np.issubdtype(rk.dtype, np.integer)
        ):
            try:
                return _read_block_fast(paths[s], rk, col_maps[s], layer=layer,
                                        as_dense=as_dense)
            except Exception:
                pass  # fall back to the robust reader on any fast-path issue
        return _read_source_block(
            paths[s], rk, out_var_names, layer=layer, as_dense=as_dense
        )

    max_workers = max(1, int(thread))
    if max_workers == 1 or len(runs) <= 1:
        for start, end, s, row_key in runs:
            write_block(start, end, _read(s, row_key))
        return

    with ThreadPoolExecutor(max_workers=max_workers) as exe:
        pending = deque()
        ri = 0
        while ri < len(runs) and len(pending) < max_workers:
            st, en, s, rk = runs[ri]
            pending.append((st, en, exe.submit(_read, s, rk)))
            ri += 1
        while pending:
            st, en, fut = pending.popleft()
            write_block(st, en, fut.result())
            if ri < len(runs):
                st2, en2, s2, rk2 = runs[ri]
                pending.append((st2, en2, exe.submit(_read, s2, rk2)))
                ri += 1



def _read_source_block(
    path: str,
    row_key: np.ndarray,
    out_var_names: List[str],
    layer: Optional[str] = None,
    as_dense: bool = False,
) -> sp.csr_matrix:
    """Read a block of rows from one source `.h5ad` as a CSR matrix.

    Parameters
    ----------
    path : str
        Source `.h5ad` file.
    row_key : np.ndarray
        Either integer positions (`_orig_obs_pos`) or obs-name labels
        (`_orig_obs_name`) selecting rows, in the desired output order.
    out_var_names : list of str
        Target variable order; columns are reordered/selected to match.
    layer : str or None, optional
        Name of the source layer to read. When `None` (default) the main `X`
        matrix is used; otherwise the named entry in `layers` is read.
    as_dense : bool, optional
        When ``True`` return a dense float32 ndarray instead of a CSR matrix.
        Default ``False``.

    Returns
    -------
    scipy.sparse.csr_matrix or numpy.ndarray
        Shape `(len(row_key), len(out_var_names))`, dtype float32. Dense when
        `as_dense` is True, otherwise CSR.
    """
    a = anndata.read_h5ad(path, backed="r")
    try:
        row_key = np.asarray(row_key)
        if np.issubdtype(row_key.dtype, np.integer):
            X_rows = _read_rows_in_order(_matrix_row_source(a, layer), row_key)
        else:
            sub = a[row_key, :]
            X_rows = sub.X if layer is None else sub.layers[layer]
        # Normalize: sparse -> CSR for fast column slicing. A dense source is
        # kept dense when a dense result is requested, avoiding a wasteful
        # dense->CSR->dense round trip; otherwise it is densified to CSR.
        if sp.issparse(X_rows):
            X_rows = X_rows.tocsr()
        elif not as_dense:
            X_rows = sp.csr_matrix(X_rows)
        else:
            X_rows = np.asarray(X_rows)

        # Map source var positions to the output var order.
        src_var = list(map(str, a.var_names))
        n_out = len(out_var_names)
        # Fast path: source var order already matches the output order
        # (the common case — `from_files` requires identical vars). Avoids a
        # full CSR column-gather copy on every block.
        if src_var == list(out_var_names):
            out = X_rows
        else:
            pos_map = {v: i for i, v in enumerate(src_var)}
            col_pos = np.array([pos_map.get(v, -1) for v in out_var_names], dtype=int)
            present = col_pos >= 0
            if present.all():
                out = X_rows[:, col_pos]
            else:
                out = sp.lil_matrix((X_rows.shape[0], n_out), dtype=np.float32)
                out[:, np.nonzero(present)[0]] = X_rows[:, col_pos[present]]
        if as_dense:
            if sp.issparse(out):
                out = out.toarray()
            return np.asarray(out, dtype=np.float32)
        if not sp.isspmatrix_csr(out):
            out = out.tocsr()
        return out.astype(np.float32, copy=False)
    finally:
        _close_backed(a)



def main():
    """CLI entry point using Python Fire."""
    import fire

    fire.core.Display = lambda lines, out: print(*lines, file=out)  # type: ignore
    fire.Fire()


if __name__ == "__main__":
    main()
