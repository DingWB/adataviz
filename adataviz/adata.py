import os
from typing import List, Dict, Optional, Sequence, Tuple, Any
import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import glob
import logging
import json

log = logging.getLogger(__name__)

class AnnDataCollection:
    """Wrapper representing a merged collection of AnnData files.

    Example:
        ds = AnnDataCollection.from_files(["a.h5ad", "b.h5ad"], out_path="merged.h5ad")
        view = ds[[True]*100, ["GeneA", "GeneB"]]
        adata = view.to_memory()
    """

    def __init__(self, adata_paths: List[str], 
                 ann: Optional[anndata.AnnData] = None, 
                 source_info: Optional[List[Dict]] = None):
        """Initialize an `AnnDataCollection` wrapper.

        Parameters
        - `adata_paths`: list of filesystem paths to the individual `.h5ad` files.
        - `ann`: optional `anndata.AnnData` that contains merged metadata
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
        self.ann = ann
        self.source_info = source_info or []

    @classmethod
    def from_files(cls, paths: Sequence[str], 
                   out_path: Optional[str] = None) -> "AnnDataCollection":
        """Create a merged AnnDataCollection from existing `.h5ad` files.
        This will merge `obs` (stacked) and `var` (union by var_names).
        `X` is not merged; the saved on-disk AnnData will contain an empty
        sparse matrix of shape (n_obs_total, n_vars_total).
        """
        # `paths` should be an iterable of paths to `.h5ad` files. When
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
        var_indexes = pd.Index([])
        per_source_info: List[Dict[str, Any]] = []
        n_obs_total = 0

        # First pass: gather var union, obs frames and per-source info.
        # Read in backed mode so we don't load `X` into memory.
        for i, p in enumerate(paths):
            a = anndata.read_h5ad(p, backed="r")
            try:
                n_obs = a.n_obs
                obs = a.obs.copy()
                # avoid reset_index which can trigger anndata's
                # ImplicitModificationWarning by directly adding a column
                # with the original obs names (as strings)
                obs["_orig_obs_name"] = a.obs_names.astype(str)
                obs["_source_idx"] = i

                obs_list.append(obs)
                var_indexes = var_indexes.union(a.var_names)

                per_source_info.append({
                    "path": os.path.abspath(p),
                    "n_obs": n_obs,
                    "var_names": list(a.var_names),
                })

                n_obs_total += n_obs
            finally:
                try:
                    a.file.close()
                except Exception:
                    try:
                        # fallback for older anndata versions
                        a._file.close()
                    except Exception:
                        pass

        # Concatenate obs
        merged_obs = pd.concat(obs_list, ignore_index=True)
        # ensure obs index is string-typed to avoid anndata implicit index conversion warning
        merged_obs.index = merged_obs.index.map(str)

        # Build merged var DataFrame (union of var names)
        merged_var = pd.DataFrame(index=var_indexes)
        # Create an empty sparse X with proper shape (we don't merge expression matrices here)
        X_empty = sp.csr_matrix((n_obs_total, len(merged_var)), dtype=np.float32)

        merged = anndata.AnnData(X_empty, obs=merged_obs, var=merged_var)
        # Note: we intentionally do not merge `obsm` arrays. If per-source
        # `obsm` data is needed, read individual files directly from
        # `individual_adata_paths` stored in `merged.uns`.

        # store metadata
        merged.uns["individual_adata_paths"] = [os.path.abspath(p) for p in paths]
        # JSON-serialize per-source info so HDF5 can store it safely
        try:
            merged.uns["_adataset_source_info"] = json.dumps(per_source_info)
        except Exception:
            # fallback: store as a list of strings
            merged.uns["_adataset_source_info"] = [str(x) for x in per_source_info]

        instance = cls(list(paths), ann=merged, source_info=per_source_info)
        # optionally save to out_path
        if out_path:
            merged.write_h5ad(out_path)

        return instance

    def __len__(self) -> int:
        return self.ann.n_obs if self.ann is not None else 0

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

    def _parse_obs_indexer(self, idx) -> np.ndarray:
        n = self.ann.n_obs
        if idx is None:
            return np.arange(n)
        if isinstance(idx, slice):
            return np.arange(n)[idx]
        if isinstance(idx, (list, np.ndarray, pd.Series)):
            arr = np.asarray(idx)
            if arr.dtype == bool:
                return np.nonzero(arr)[0]
            # If integer dtype, treat as indices. Otherwise treat as labels.
            if np.issubdtype(arr.dtype, np.integer):
                return arr.astype(int)
            # treat as obs names (strings or object)
            return np.array(self.ann.obs_names.get_indexer(arr.astype(str)), dtype=int)
        raise IndexError("Unsupported obs indexer")

    def _parse_var_indexer(self, idx) -> List[str]:
        if idx is None:
            return list(self.ann.var_names)
        if isinstance(idx, slice):
            return list(self.ann.var_names[idx])
        if isinstance(idx, (list, np.ndarray, pd.Series)):
            arr = np.asarray(idx)
            if arr.dtype == bool:
                return list(np.asarray(self.ann.var_names)[arr])
            # If integer dtype, interpret as indices into var_names
            if np.issubdtype(arr.dtype, np.integer):
                return list(np.asarray(self.ann.var_names)[arr.astype(int)])
            # Otherwise treat as labels (string/object/unicode)
            return list(arr.astype(str))
        if isinstance(idx, str):
            return [idx]
        raise IndexError("Unsupported var indexer")

class AnnDataView:
    """Helper representing a (cells, genes) subset of an `AnnDataCollection`.

    Call `to_memory()` to read and assemble the actual `AnnData` with `X`.
    """

    def __init__(self, dataset: AnnDataCollection, obs_idx: np.ndarray, var_names: List[str]):
        self.dataset = dataset
        self.obs_idx = np.asarray(obs_idx, dtype=int)
        self.var_names = list(var_names)

    def to_memory(self) -> anndata.AnnData:
        """Load selected `X` chunks from underlying files and assemble AnnData.

        Returns a new `anndata.AnnData` with concatenated `X` for the selected
        cells and genes.
        """
        ds = self.dataset
        merged = ds.ann

        # Build target obs and var
        sel_obs = merged.obs.iloc[self.obs_idx].copy()
        sel_var = merged.var.loc[self.var_names].copy()

        n_rows = len(self.obs_idx)
        n_cols = len(self.var_names)

        # Prepare a sparse matrix to fill
        X_global = sp.lil_matrix((n_rows, n_cols), dtype=np.float32)

        # Map merged obs rows to sources
        src_indices = merged.obs.iloc[self.obs_idx]["_source_idx"].values
        # For each source, find positions and corresponding original obs names
        for src_idx, path in enumerate(ds.adata_paths):
            mask = src_indices == src_idx
            if not np.any(mask):
                continue

            # global row positions that come from this source
            global_rows = np.nonzero(mask)[0]
            # original obs names stored in `_orig_obs_name` column
            orig_obs_names = sel_obs.iloc[mask]["_orig_obs_name"].values

            a = anndata.read_h5ad(path, backed="r")
            try:
                # map orig_obs_names to positions in source adata
                src_row_pos = np.asarray(a.obs_names.get_indexer(orig_obs_names), dtype=int)
                # ensure all requested obs were found in the source
                if (src_row_pos < 0).any():
                    missing = orig_obs_names[src_row_pos < 0]
                    raise KeyError(f"Some requested obs names not found in {path}: {missing}")

                # find intersection of requested genes with this source's genes
                src_gene_mask = np.isin(self.var_names, a.var_names)
                if not np.any(src_gene_mask):
                    continue

                genes_in_source = [g for g, m in zip(self.var_names, src_gene_mask) if m]
                src_col_pos_in_global = np.nonzero(src_gene_mask)[0]

                # map gene names to source column indices
                src_gene_pos = np.asarray(a.var_names.get_indexer(genes_in_source), dtype=int)
                # sanity check: ensure genes_in_source were located
                if (src_gene_pos < 0).any():
                    missing_g = np.asarray(genes_in_source)[src_gene_pos < 0]
                    raise KeyError(f"Some requested genes not found in {path}: {missing_g}")

                # read X slice from source. Some anndata/h5ad backends may not
                # support simultaneous row+col indexing in `backed='r'` mode,
                # so try direct slicing and fall back to row-first then
                # column-subsetting if needed.
                try:
                    sub = a[src_row_pos, src_gene_pos]
                    X_sub = sub.X
                except Exception:
                    # fallback: select rows first, then subset columns
                    rows_ann = a[src_row_pos]
                    X_rows = rows_ann.X
                    if sp.issparse(X_rows):
                        X_rows = X_rows.tocsr()
                    else:
                        X_rows = sp.csr_matrix(X_rows)
                    X_sub = X_rows[:, src_gene_pos]

                if sp.issparse(X_sub):
                    X_sub = X_sub.tocsr()
                else:
                    X_sub = sp.csr_matrix(X_sub)

                # place into global
                for i_local, i_global in enumerate(global_rows):
                    row = X_sub[i_local]
                    if sp.issparse(row):
                        row = row.toarray().ravel()
                    else:
                        row = np.asarray(row).ravel()
                    if row.size != len(src_col_pos_in_global):
                        # unexpected shape
                        raise RuntimeError("Shape mismatch when assembling X")
                    X_global[i_global, src_col_pos_in_global] = row
            finally:
                try:
                    a.file.close()
                except Exception:
                    try:
                        a._file.close()
                    except Exception:
                        pass

        X_final = X_global.tocsr()
        out = anndata.AnnData(X_final, obs=sel_obs.reset_index(drop=True), var=sel_var)
        return out

def _compare_annodata_metadata(a1: anndata.AnnData, a2: anndata.AnnData) -> Tuple[bool, str]:
    """Compare metadata (obs and var) of two AnnData objects.

    Returns (equal, message). Does not compare `X`.
    """
    import pandas.testing as pdt

    try:
        if a1.n_obs != a2.n_obs:
            return False, f"n_obs differ: {a1.n_obs} != {a2.n_obs}"
        if a1.n_vars != a2.n_vars:
            return False, f"n_vars differ: {a1.n_vars} != {a2.n_vars}"

        # compare names order
        if list(a1.obs_names) != list(a2.obs_names):
            return False, "obs_names differ"
        if list(a1.var_names) != list(a2.var_names):
            return False, "var_names differ"

        # compare obs DataFrame
        try:
            pdt.assert_frame_equal(a1.obs, a2.obs, check_dtype=False, check_like=False)
        except AssertionError as e:
            return False, f"obs DataFrame mismatch: {e}"

        # compare var DataFrame
        try:
            pdt.assert_frame_equal(a1.var, a2.var, check_dtype=False, check_like=False)
        except AssertionError as e:
            return False, f"var DataFrame mismatch: {e}"

        return True, "metadata identical"
    except Exception as e:
        return False, f"comparison failed: {e}"

def _compare_adata_overlap(ds: "AnnDataCollection", reference_path: str, tol: float = 1e-6) -> Tuple[bool, str]:
    """Compare expression values for overlapping obs and var between a
    merged `AnnDataCollection` and a reference `.h5ad` file.
    - `ds` is an `AnnDataCollection` instance (metadata merged, X assembled on demand).
    - `reference_path` is the path to the merged reference `.h5ad` file.
    Returns (equal, message).
    """
    ref = anndata.read_h5ad(reference_path, backed="r")
    try:
        # compute intersections while preserving order from the reference
        obs_common = np.intersect1d(np.asarray(ds.ann.obs_names, dtype=str), np.asarray(ref.obs_names, dtype=str), assume_unique=False)
        var_common = np.intersect1d(np.asarray(ds.ann.var_names, dtype=str), np.asarray(ref.var_names, dtype=str), assume_unique=False)

        if len(obs_common) == 0 or len(var_common) == 0:
            return False, "No overlapping obs or vars to compare"

        # assemble dataset X for the overlapping obs/vars
        subset = ds[list(obs_common), list(var_common)]
        assembled = subset.to_memory()

        # read the corresponding slice from the reference (backed)
        ref_obs_idx = np.asarray(ref.obs_names.get_indexer(obs_common), dtype=int)
        ref_var_idx = np.asarray(ref.var_names.get_indexer(var_common), dtype=int)
        try:
            ref_sub = ref[ref_obs_idx, ref_var_idx]
            X_ref = ref_sub.X
        except Exception:
            tmp = ref[ref_obs_idx]
            X_tmp = tmp.X
            if sp.issparse(X_tmp):
                X_tmp = X_tmp.tocsr()
            else:
                X_tmp = sp.csr_matrix(X_tmp)
            X_ref = X_tmp[:, ref_var_idx]

        X_a = assembled.X
        X_b = X_ref

        # convert to csr
        if not sp.issparse(X_a):
            X_a = sp.csr_matrix(X_a)
        else:
            X_a = X_a.tocsr()
        if not sp.issparse(X_b):
            X_b = sp.csr_matrix(X_b)
        else:
            X_b = X_b.tocsr()

        if X_a.shape != X_b.shape:
            return False, f"shape mismatch for overlap: {X_a.shape} vs {X_b.shape}"

        diff = X_a - X_b
        # compute maximum absolute difference
        if sp.issparse(diff):
            max_abs = float(np.max(np.abs(diff.data))) if diff.nnz > 0 else 0.0
        else:
            max_abs = float(np.max(np.abs(diff)))

        if np.isfinite(max_abs) and max_abs <= tol:
            return True, f"X equal within tol={tol} (max_abs={max_abs})"
        else:
            return False, f"X differs: max_abs={max_abs} > tol={tol}"
    finally:
        try:
            ref.file.close()
        except Exception:
            try:
                ref._file.close()
            except Exception:
                pass

def test_merge_and_subset(src_path: str, out_dir: str, 
                          n_obs_subset: int = 10, 
                          n_var_subset: int = 20,
                          reference_path: Optional[str] = None):
    """Simple test utility: merge up to `max_files` .h5ad files from `src_dir`,
    write merged metadata to `out_dir`, then test `AnnDataView.to_memory()` by
    selecting the first `n_obs_subset` observations and first `n_var_subset`
    variables and assembling the expression matrix.

    Prints status and shapes; raises on failure.
    """
    files = sorted(glob.glob(os.path.expanduser(src_path)))

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "merged_test_metadata.h5ad")

    log.info("Merging %d files", len(files))
    ds = AnnDataCollection.from_files(files, out_path=out_path)
    log.info("Merged metadata saved to %s", out_path)

    # If a reference merged file is provided, compare metadata
    if reference_path:
        ref = anndata.read_h5ad(reference_path, backed="r")
        try:
            equal, msg = _compare_annodata_metadata(ds.ann, ref)
            if equal:
                log.info("Merged dataset metadata is IDENTICAL to reference: %s", reference_path)
            else:
                log.warning("Merged dataset metadata differs from reference: %s", msg)
        finally:
            try:
                ref.file.close()
            except Exception:
                try:
                    ref._file.close()
                except Exception:
                    pass

    # choose the first n_obs_subset rows and first n_var_subset genes
    n_obs = min(n_obs_subset, ds.ann.n_obs)
    if n_obs <= 0:
        raise RuntimeError("No observations in merged AnnData")
    var_names = list(ds.ann.var_names[:n_var_subset])
    obs_indexer = list(range(n_obs))

    subset = ds[obs_indexer, var_names]
    log.info("Assembling subset to memory (obs=%d, vars=%d)", n_obs, len(var_names))
    assembled = subset.to_memory()
    log.info("Assembled AnnData: obs=%d, var=%d, X shape=%s", assembled.n_obs, assembled.n_vars, getattr(assembled.X, 'shape', None))
    return assembled

def run_test():
    assembled = test_merge_and_subset("/anvil/projects/x-mcb130189/Wubin/CB/adata/100kb/*-CGN.h5ad", 
                               "/anvil/projects/x-mcb130189/Wubin/test", 
                               n_obs_subset=10, n_var_subset=20,
                               reference_path="/anvil/projects/x-mcb130189/Wubin/CB/adata/CB_BICAN.100kb-CGN.h5ad")
    print("Test completed:", assembled.shape if hasattr(assembled, 'shape') else (assembled.n_obs, assembled.n_vars))

def main():
	import fire
	fire.core.Display = lambda lines, out: print(*lines, file=out) # type: ignore
	fire.Fire()

if __name__ == "__main__":
	main()
