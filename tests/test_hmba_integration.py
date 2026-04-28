"""Integration tests against a real HMBA AnnData file.

These tests are skipped automatically when the data file or required
optional dependencies are not available. To run them explicitly::

    /path/to/m3c/bin/python -m pytest tests/test_hmba_integration.py -v
"""

import os
import tempfile

import matplotlib

matplotlib.use("Agg")  # noqa: E402  (must be set before pyplot import)

import anndata
import matplotlib.pyplot as plt
import pandas as pd
import pytest

HMBA_PATH = os.path.expanduser(
    "~/Projects/LabData/BasalGanglia/10xMultiome/HMBA.Group.downsample_1500.h5ad"
)

pytestmark = pytest.mark.skipif(
    not os.path.exists(HMBA_PATH),
    reason=f"HMBA test fixture not found at {HMBA_PATH}",
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def hmba_adata_backed():
    """Backed read of the HMBA AnnData. Closed at module teardown."""
    adata = anndata.read_h5ad(HMBA_PATH, backed="r")
    yield adata
    try:
        adata.file.close()
    except Exception:
        pass


@pytest.fixture(scope="module")
def hmba_obs():
    """Cached copy of ``adata.obs`` for tests that don't need ``X``."""
    adata = anndata.read_h5ad(HMBA_PATH, backed="r")
    obs = adata.obs.copy()
    adata.file.close()
    return obs


@pytest.fixture
def tmp_outdir(tmp_path):
    return str(tmp_path)


# ---------------------------------------------------------------------------
# AnnDataCollection
# ---------------------------------------------------------------------------


class TestAnnDataCollection:
    def test_from_files_single_path(self, tmp_outdir, hmba_adata_backed):
        from adataviz.adata import AnnDataCollection, is_anndatacollection

        out = os.path.join(tmp_outdir, "merged.h5ad")
        ds = AnnDataCollection.from_files([HMBA_PATH], out_path=out)
        assert os.path.exists(out)
        assert ds.adata.n_obs == hmba_adata_backed.n_obs
        assert ds.adata.n_vars == hmba_adata_backed.n_vars
        assert is_anndatacollection(out) is True

    def test_subset_cells_matches_reference(self, tmp_outdir, hmba_adata_backed):
        from adataviz.adata import AnnDataCollection

        out = os.path.join(tmp_outdir, "merged.h5ad")
        ds = AnnDataCollection.from_files([HMBA_PATH], out_path=out)

        rng = pd.Series(hmba_adata_backed.obs_names)
        use_cells = rng.sample(200, random_state=0).tolist()
        use_genes = (
            pd.Series(hmba_adata_backed.var_names).sample(20, random_state=0).tolist()
        )

        view = ds[use_cells, use_genes]
        sub = view.to_memory(thread=2)
        assert sub.n_obs == 200
        assert sub.n_vars == 20
        # Order of var_names is preserved as requested.
        assert sub.var_names.tolist() == use_genes


# ---------------------------------------------------------------------------
# adataviz.tl (tools)
# ---------------------------------------------------------------------------


class TestTools:
    def test_get_obs(self, tmp_outdir):
        from adataviz import tl

        outfile = os.path.join(tmp_outdir, "obs.tsv")
        path = tl.get_obs(HMBA_PATH, add_coord=True, outfile=outfile)
        assert path == outfile
        df = pd.read_csv(outfile, sep="\t", index_col=0)
        assert "Subclass" in df.columns
        assert "Group" in df.columns
        # add_coord=True should have added umap_0 / umap_1 columns
        assert any(c.startswith("umap") for c in df.columns)
        assert df.shape[0] > 0

    def test_downsample_adata(self, tmp_outdir, hmba_obs):
        from adataviz import tl

        outfile = os.path.join(tmp_outdir, "downsampled.h5ad")
        tl.downsample_adata(
            HMBA_PATH, groupby="Subclass", outfile=outfile, downsample=20
        )
        assert os.path.exists(outfile)
        ad = anndata.read_h5ad(outfile, backed="r")
        try:
            # Each Subclass should now have at most 20 cells.
            counts = ad.obs["Subclass"].value_counts()
            assert (counts <= 20).all()
            assert ad.n_obs <= hmba_obs.shape[0]
        finally:
            ad.file.close()

    def test_composition(self, tmp_outdir, hmba_obs):
        from adataviz import tl

        cwd = os.getcwd()
        os.chdir(tmp_outdir)
        try:
            obs = hmba_obs.copy()
            tl.composition(
                obs=obs,
                groupby="Subclass",
                composition_col="Class",
                outname="composition.xlsx",
            )
            outfile = os.path.join(tmp_outdir, "composition.xlsx")
            assert os.path.exists(outfile)
            df = pd.read_excel(outfile, sheet_name="All", index_col=0)
            # rows should be Subclass categories present in obs
            assert df.shape[0] > 0
        finally:
            os.chdir(cwd)


# ---------------------------------------------------------------------------
# adataviz.pl (plotting)
# ---------------------------------------------------------------------------


class TestPlotting:
    def test_plot_categorical_subclass(self, tmp_outdir, hmba_adata_backed):
        from adataviz import pl

        # Use a small in-memory subset so the plot is fast.
        sub = hmba_adata_backed[:2000, :].to_memory()
        save = os.path.join(tmp_outdir, "umap_subclass.png")
        pl.plot_categorical(
            sub,
            basis="umap",
            groupby="Subclass",
            save=save,
            show=False,
            figsize=(4, 3),
        )
        plt.close("all")
        assert os.path.exists(save)
        assert os.path.getsize(save) > 0

    def test_plot_categorical_group(self, tmp_outdir, hmba_adata_backed):
        from adataviz import pl

        sub = hmba_adata_backed[:2000, :].to_memory()
        save = os.path.join(tmp_outdir, "umap_group.png")
        pl.plot_categorical(
            sub,
            basis="umap",
            groupby="Group",
            save=save,
            show=False,
            figsize=(4, 3),
            coding=False,
        )
        plt.close("all")
        assert os.path.exists(save)
