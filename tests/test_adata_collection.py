"""Integration tests for ``adataviz.adata.AnnDataCollection``.

These tests rely on local on-disk data and are skipped by default. Run with::

    ADATAVIZ_RUN_INTEGRATION=1 pytest tests/test_adata_collection.py
"""

import os
import time

import anndata
import pandas as pd
import pytest

RUN_INTEGRATION = os.environ.get("ADATAVIZ_RUN_INTEGRATION") == "1"


@pytest.mark.skipif(
    not RUN_INTEGRATION, reason="set ADATAVIZ_RUN_INTEGRATION=1 to enable"
)
def test_anndatacollection_subset():
    """Compare AnnDataCollection subset against direct backed AnnData read."""
    from loguru import logger
    from adataviz.adata import AnnDataCollection

    os.chdir(
        os.path.expanduser("~/Projects/mouse_dev/testAnnDataCollection")
    )
    adata_path = "~/Projects/mouse_dev/adata/100kb/*-CGN.h5ad"
    reference_path = (
        "~/Projects/mouse_dev/adata/mouse_dev.100kb-CGN.h5ad"
    )
    metadata_path = os.path.expanduser(
        "~/Projects/mouse_dev/metadata/metadata.filtered.tsv.gz"
    )
    ds = AnnDataCollection.from_files(
        adata_path,
        out_path="mouse_dev.100kb-CGN.h5ad",
        metadata_path=metadata_path,
    )

    ref = anndata.read_h5ad(reference_path, backed="r")
    use_cells = ref.obs.sample(5000).index.tolist()
    start = time.perf_counter()
    ref_data = ref[use_cells, :].to_df()
    logger.info(f"Reference read time for 5000 cells: {time.perf_counter() - start:.2f}s")
    ref.file.close()

    start = time.perf_counter()
    data = ds[use_cells, :].to_memory(thread=10).to_df()
    logger.info(
        f"AnnDataCollection read time for 5000 cells: {time.perf_counter() - start:.2f}s"
    )
    ref_data.index.name = "cell"
    pd.testing.assert_frame_equal(ref_data, data, check_dtype=True, rtol=1e-6, atol=0)

    ref = anndata.read_h5ad(reference_path, backed="r")
    use_vars = ref.var.sample(50).index.tolist()
    ref_data = ref[:, use_vars].to_df()
    ref.file.close()

    data = ds[:, use_vars].to_memory(thread=10).to_df()
    ref_data.index.name = "cell"
    data = data.loc[ref_data.index.tolist()]
    pd.testing.assert_frame_equal(ref_data, data, check_dtype=True, rtol=1e-6, atol=0)
