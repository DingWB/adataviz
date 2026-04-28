"""Unit tests for the cell-statistics plot functions in :mod:`adataviz.pl`."""

import os

import matplotlib

matplotlib.use("Agg")  # noqa: E402

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from adataviz import pl


@pytest.fixture
def cell_obs():
    """Synthetic cell metadata with two categorical axes and a marker set."""
    rng = np.random.default_rng(0)
    n = 600
    cell_types = np.array(["TypeA", "TypeB", "TypeC", "TypeD"])
    regions = np.array(["R1", "R2", "R3"])
    donors = np.array(["donor1", "donor2", "donor3"])
    df = pd.DataFrame(
        {
            "cell_type": rng.choice(cell_types, size=n, p=[0.4, 0.3, 0.2, 0.1]),
            "region": rng.choice(regions, size=n),
            "donor": rng.choice(donors, size=n),
            "marker": rng.choice(
                [f"gene_{i}" for i in range(15)], size=n, replace=True
            ),
        },
        index=[f"cell_{i}" for i in range(n)],
    )
    return df


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# Individual plot functions
# ---------------------------------------------------------------------------


class TestRose:
    def test_basic(self, cell_obs):
        ax = pl.rose_plot(cell_obs, groupby="cell_type")
        assert ax.name == "polar"

    def test_split(self, cell_obs):
        ax = pl.rose_plot(cell_obs, groupby="cell_type", split_by="region")
        # Each split adds patches; should have at least 4 wedges.
        assert len(ax.patches) >= 4


class TestRing:
    def test_returns_axes(self, cell_obs):
        ax = pl.ring_plot(cell_obs, groupby="cell_type")
        assert hasattr(ax, "figure")
        # Donut wedges = number of categories
        from matplotlib.patches import Wedge

        wedges = [p for p in ax.patches if isinstance(p, Wedge)]
        assert len(wedges) == cell_obs["cell_type"].nunique()

    def test_save(self, cell_obs, tmp_path):
        out = tmp_path / "ring.png"
        pl.ring_plot(cell_obs, groupby="cell_type", save=str(out))
        assert out.exists() and out.stat().st_size > 0


class TestPie:
    def test_simple(self, cell_obs):
        ax = pl.pie_plot(cell_obs, groupby="cell_type")
        assert hasattr(ax, "figure")

    def test_facets(self, cell_obs):
        axes = pl.pie_plot(cell_obs, groupby="cell_type", split_by="region")
        # at least one axes per region
        assert len(axes) >= cell_obs["region"].nunique()


class TestArea:
    def test_normalized(self, cell_obs):
        ax = pl.area_plot(cell_obs, groupby="cell_type", split_by="region")
        # Each x position should sum to 1 in stackplot
        assert ax.get_ylim()[1] == pytest.approx(1)

    def test_unnormalized(self, cell_obs):
        ax = pl.area_plot(
            cell_obs, groupby="cell_type", split_by="region", normalize=False
        )
        assert ax.get_ylabel() == "Count"


class TestDot:
    def test_basic(self, cell_obs):
        ax = pl.dot_plot(cell_obs, groupby="cell_type", split_by="region")
        # one PathCollection from scatter
        assert len(ax.collections) >= 1
        n_dots = ax.collections[0].get_offsets().shape[0]
        assert n_dots == cell_obs["cell_type"].nunique() * cell_obs["region"].nunique()


class TestTrend:
    def test_basic(self, cell_obs):
        ax = pl.trend_plot(cell_obs, groupby="cell_type", split_by="region")
        assert len(ax.lines) == cell_obs["cell_type"].nunique()


class TestSankey:
    def test_basic(self, cell_obs):
        ax = pl.sankey_plot(cell_obs, left="cell_type", right="region")
        # Nodes (rectangles) + flows (PathPatch); at least n_left + n_right rects
        from matplotlib.patches import Rectangle

        rects = [p for p in ax.patches if isinstance(p, Rectangle)]
        assert (
            len(rects)
            >= cell_obs["cell_type"].nunique() + cell_obs["region"].nunique()
        )


class TestChord:
    def test_basic(self, cell_obs):
        # Reduce the number of categories to keep the chord readable.
        sub = cell_obs.copy()
        fig = pl.chord_plot(sub, left="cell_type", right="region")
        # accept either a Figure or pycirclize-returned Figure
        assert fig is not None


class TestVenn:
    def test_two_sets(self, cell_obs):
        sub = cell_obs[cell_obs["region"].isin(["R1", "R2"])].copy()
        ax = pl.venn_plot(sub, groupby="region", set_by="marker")
        assert hasattr(ax, "figure")

    def test_three_sets(self, cell_obs):
        ax = pl.venn_plot(cell_obs, groupby="region", set_by="marker")
        assert hasattr(ax, "figure")

    def test_too_many_sets(self, cell_obs):
        with pytest.raises(ValueError):
            pl.venn_plot(cell_obs, groupby="cell_type", set_by="marker")


class TestUpset:
    def test_basic(self, cell_obs):
        fig = pl.upset_plot(cell_obs, groupby="cell_type", set_by="marker")
        assert fig is not None
