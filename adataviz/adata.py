import os
from typing import List, Dict, Optional, Sequence, Tuple, Any
from concurrent.futures import ThreadPoolExecutor, as_completed
import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import glob
import json
from loguru import logger as logger
# logger.remove()
# logger.add(sys.stderr, level="DEBUG")


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
    ) -> "AnnDataCollection":
        """Create a merged AnnDataCollection from existing `.h5ad` files.
        This will merge `obs` (stacked) and `var` (union by var_names).
        `X` is not merged; the saved on-disk AnnData will contain an empty
        sparse matrix of shape (n_obs_total, n_vars_total).
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
                try:
                    a.file.close()
                except Exception:
                    try:
                        a._file.close()
                    except Exception:
                        pass

        # First pass: gather obs frames and per-source info.
        # Read in backed mode so we don't load X into memory.
        # Parallelize across source files for faster metadata gathering.
        n_workers = min(8, max(1, len(paths)))
        results_by_idx: Dict[int, Any] = {}
        with ThreadPoolExecutor(max_workers=n_workers) as exe:
            for i, obs, info, src_vars in exe.map(_read_source_meta, enumerate(paths)):
                results_by_idx[i] = (obs, info, src_vars)
        # Process results in order to preserve deterministic obs ordering
        for i in range(len(paths)):
            obs, info, src_vars = results_by_idx[i]
            obs_list.append(obs)
            per_source_info.append(info)
            if var_names is None:
                var_names = src_vars
            else:
                assert var_names == src_vars, (
                    f"Var names in file {paths[i]} differ from previous files; this may lead to incorrect merging."
                )

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
            meta = pd.read_csv(metadata_path, index_col=0, sep=sep)
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
        try:
            ad.file.close()
        except Exception:
            try:
                ad._file.close()
            except Exception:
                pass


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

    def to_memory(self, thread=8) -> anndata.AnnData:
        """Load selected `X` chunks from underlying files and assemble AnnData.
        Returns a new `anndata.AnnData` with concatenated `X` for the selected
        cells and genes.
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

        # Check for integer position column (faster HDF5 access)
        has_pos = "_orig_obs_pos" in sel_obs.columns
        if has_pos:
            sel_orig_pos = sel_obs["_orig_obs_pos"].values.astype(int)
        sel_orig_obs = sel_obs["_orig_obs_name"].to_numpy()

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
        src_col_pos_in_output = np.nonzero(present_mask)[0]
        # Positions in each *source* file's var axis
        src_gene_pos = src_gene_pos_all[present_mask]

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

        def _process_source(item):
            """Read and extract X data from a single source h5ad file."""
            src_idx, path, global_rows, row_key = item
            a = anndata.read_h5ad(path, backed="r")
            try:
                # Read selected rows from disk
                if np.issubdtype(row_key.dtype, np.integer):
                    # Sort positions for sequential HDF5 reads, then restore original order
                    sort_order = np.argsort(row_key)
                    sorted_pos = row_key[sort_order]
                    unsort_order = np.argsort(sort_order)
                    X_rows = a.X[sorted_pos, :]
                    # Restore original row order
                    if sp.issparse(X_rows):
                        X_rows = X_rows.tocsr()[unsort_order, :]
                    else:
                        X_rows = X_rows[unsort_order, :]
                else:
                    X_rows = a[row_key, :].X

                # Ensure CSR for fast column slicing
                if not sp.issparse(X_rows):
                    X_rows = sp.csr_matrix(X_rows)
                elif not sp.isspmatrix_csr(X_rows):
                    X_rows = X_rows.tocsr()

                # Extract only the needed columns
                X_sub = X_rows[:, src_gene_pos]
                X_coo = sp.coo_matrix(X_sub)

                if X_coo.nnz == 0:
                    return None

                rows_mapped = global_rows[X_coo.row].astype(np.int32)
                cols_mapped = src_col_pos_in_output[X_coo.col].astype(np.int32)
                data_arr = X_coo.data.astype(np.float32)

                return rows_mapped, cols_mapped, data_arr
            finally:
                try:
                    a.file.close()
                except Exception:
                    try:
                        a._file.close()
                    except Exception:
                        pass

        # Process only sources that have selected cells, in parallel
        rows_all: List[np.ndarray] = []
        cols_all: List[np.ndarray] = []
        data_all: List[np.ndarray] = []

        max_workers = min(thread, max(1, len(work_items)))
        if max_workers == 0:
            pass  # no work to do
        elif max_workers == 1:
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


# Backward-compatible alias for the previous (mis-spelled) public name.
is_annadatacollection = is_anndatacollection


def main():
    """CLI entry point using Python Fire."""
    import fire

    fire.core.Display = lambda lines, out: print(*lines, file=out)  # type: ignore
    fire.Fire()


if __name__ == "__main__":
    main()
