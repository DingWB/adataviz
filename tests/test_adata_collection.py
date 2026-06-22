"""Integration tests for ``adataviz.adata.AnnDataCollection``.

These tests rely on local on-disk data. Run with::

    pytest tests/test_adata_collection.py
"""

import os
import threading
import time
import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _rss_bytes():
    """Current resident set size (RSS) of this process, in bytes (Linux)."""
    with open("/proc/self/statm") as f:
        rss_pages = int(f.read().split()[1])
    return rss_pages * os.sysconf("SC_PAGE_SIZE")


class peak_memory:
    """Context manager that samples RSS in a background thread and reports
    the peak memory observed within the ``with`` block.

    Example
    -------
    >>> with peak_memory() as pm:
    ...     do_work()
    >>> print(pm.peak_mb)
    """

    def __init__(self, interval=0.05):
        self.interval = interval
        self.start_rss = 0
        self.peak = 0
        self._stop = threading.Event()
        self._thread = None

    def _run(self):
        while not self._stop.is_set():
            try:
                rss = _rss_bytes()
                if rss > self.peak:
                    self.peak = rss
            except Exception:
                pass
            self._stop.wait(self.interval)

    def __enter__(self):
        self.start_rss = self.peak = _rss_bytes()
        self._stop.clear()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
        rss = _rss_bytes()  # final sample
        if rss > self.peak:
            self.peak = rss
        return False

    @property
    def peak_mb(self):
        """Peak RSS observed in the block, in MiB."""
        return self.peak / (1024 * 1024)

    @property
    def increase_mb(self):
        """Increase in RSS from block entry to peak, in MiB."""
        return (self.peak - self.start_rss) / (1024 * 1024)


# Local on-disk data locations used by the integration tests.
WORKDIR = os.path.expanduser("~/Projects/mouse_dev/testAnnDataCollection")
ADATA_GLOB = "~/Projects/mouse_dev/adata/100kb/*-CGN.h5ad"
REFERENCE_PATH = os.path.expanduser(
    "~/Projects/mouse_dev/adata/mouse_dev.100kb-CGN.h5ad"
)
METADATA_PATH = os.path.expanduser(
    "~/Projects/mouse_dev/metadata/metadata.filtered.tsv.gz"
)


def _build_collection(out_path="mouse_dev.100kb-CGN.h5ads"):
    """Build the merged AnnDataCollection from the local source files."""
    from adataviz.adata import AnnDataCollection

    os.chdir(WORKDIR)
    return AnnDataCollection.from_files(
        ADATA_GLOB,
        out_path=out_path,
        metadata_path=METADATA_PATH,
    )


def _as_csr(X):
    """Return X as a float32 CSR matrix."""
    if sp.issparse(X):
        return X.tocsr().astype(np.float32, copy=False)
    return sp.csr_matrix(np.asarray(X, dtype=np.float32))


def _assert_x_matches_reference(result, reference_path, rtol=1e-6, atol=1e-5):
    """Assert ``result.X`` equals the reference X for the same cells/genes.

    The reference is sliced to ``result``'s obs/var order so the two matrices
    are aligned element-wise before comparison.
    """
    cells = result.obs_names.tolist()
    genes = result.var_names.tolist()

    # h5py (backed mode) allows only ONE fancy index at a time, so we cannot
    # do ref[cells, genes] in one shot. Subset rows while backed (single fancy
    # index), materialize, then reorder columns in memory.
    ref = anndata.read_h5ad(reference_path, backed="r")
    try:
        ref_sub = ref[cells, :].to_memory()
    finally:
        ref.file.close()
    ref_sub = ref_sub[:, genes]

    Xr = _as_csr(ref_sub.X)
    Xf = _as_csr(result.X)
    assert Xr.shape == Xf.shape, f"shape mismatch: {Xr.shape} vs {Xf.shape}"

    # Sparse element-wise comparison (memory friendly, no dense materialization).
    delta = (Xr - Xf).tocoo()
    if delta.nnz:
        ref_vals = np.asarray(Xr[delta.row, delta.col]).ravel()
        max_abs = np.abs(delta.data).max()
        max_allowed = atol + rtol * np.abs(ref_vals).max()
        assert max_abs <= max_allowed, (
            f"X differs from reference: max abs diff {max_abs} > {max_allowed}"
        )


def test_anndatacollection_subset():
    """Compare AnnDataCollection subsets and the full X against the reference."""
    from loguru import logger

    with peak_memory() as pm:
        ds = _build_collection()

        # --- cell subset: 5000 random cells, all genes ---
        ref = anndata.read_h5ad(REFERENCE_PATH, backed="r")
        use_cells = ref.obs.sample(5000).index.tolist()
        start = time.perf_counter()
        ref_data = ref[use_cells, :].to_df()
        logger.info(
            f"Reference read time for 5000 cells: {time.perf_counter() - start:.2f}s"
        )
        ref.file.close()

        start = time.perf_counter()
        data = ds[use_cells, :].to_memory(thread=10).to_df()
        logger.info(
            f"AnnDataCollection read time for 5000 cells: {time.perf_counter() - start:.2f}s"
        )
        ref_data.index.name = "cell"
        pd.testing.assert_frame_equal(ref_data, data, check_dtype=True, rtol=1e-6, atol=0)

        # --- gene subset: 50 random genes, all cells ---
        ref = anndata.read_h5ad(REFERENCE_PATH, backed="r")
        use_vars = ref.var.sample(50).index.tolist()
        ref_data = ref[:, use_vars].to_df()
        ref.file.close()

        data = ds[:, use_vars].to_memory(thread=10).to_df()
        ref_data.index.name = "cell"
        data = data.loc[ref_data.index.tolist()]
        pd.testing.assert_frame_equal(ref_data, data, check_dtype=True, rtol=1e-6, atol=0)

        # --- full X: every cell and gene must match the reference exactly ---
        start = time.perf_counter()
        full = ds[:, :].to_memory(thread=10)
        logger.info(
            f"AnnDataCollection full to_memory time: {time.perf_counter() - start:.2f}s"
        )
        _assert_x_matches_reference(full, REFERENCE_PATH)

    logger.info(
        f"[test_anndatacollection_subset] peak memory: {pm.peak_mb:.1f} MiB "
        f"(+{pm.increase_mb:.1f} MiB over baseline)"
    )


def test_anndatacollection_write_h5ad():
    """AnnDataCollection.write_h5ad must produce an h5ad whose X matches the reference."""
    from loguru import logger

    with peak_memory() as pm:
        ds = _build_collection()

        out_path = "mouse_dev.100kb-CGN.h5ad"
        start = time.perf_counter()
        ds.write_h5ad(out_path, thread=10)
        logger.info(
            f"AnnDataCollection.write_h5ad time: {time.perf_counter() - start:.2f}s"
        )

        written = anndata.read_h5ad(out_path)

        # obs/var alignment with the collection metadata.
        assert list(written.var_names) == list(ds.adata.var_names)
        assert list(written.obs_names) == list(ds.adata.obs_names)
        assert written.shape == (ds.adata.n_obs, ds.adata.n_vars)

        # Written X must match the reference for the same cells/genes.
        _assert_x_matches_reference(written, REFERENCE_PATH)

    logger.info(
        f"[test_anndatacollection_write_h5ad] peak memory: {pm.peak_mb:.1f} MiB "
        f"(+{pm.increase_mb:.1f} MiB over baseline)"
    )


# Run test:
# python -m pytest tests/test_adata_collection.py -v -s
