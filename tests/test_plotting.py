import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from adataviz.utils import categorical_scatter, continuous_scatter


@pytest.fixture(autouse=True)
def close_figures():
    """Close all matplotlib figures after each test to prevent state leakage."""
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# categorical_scatter (utils)
# ---------------------------------------------------------------------------


class TestCategoricalScatter:
    def test_returns_dataframe(self, scatter_df):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df, ax=ax, hue="cell_type", show_legend=False
        )
        assert isinstance(result, pd.DataFrame)

    def test_with_adata(self, adata_with_colors):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            adata_with_colors, ax=ax, basis="umap", hue="cell_type", show_legend=False
        )
        assert isinstance(result, pd.DataFrame)
        assert "x" in result.columns
        assert "y" in result.columns

    def test_result_has_hue_column(self, scatter_df):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df, ax=ax, hue="cell_type", show_legend=False
        )
        assert "hue" in result.columns

    def test_dict_palette(self, scatter_df):
        palette = {
            "TypeA": "#e41a1c",
            "TypeB": "#377eb8",
            "TypeC": "#4daf4a",
            "TypeD": "#984ea3",
            "TypeE": "#ff7f00",
        }
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df, ax=ax, hue="cell_type", palette=palette, show_legend=False
        )
        assert isinstance(result, pd.DataFrame)

    def test_str_palette(self, scatter_df):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df, ax=ax, hue="cell_type", palette="tab10", show_legend=False
        )
        assert isinstance(result, pd.DataFrame)

    def test_with_text_annotation(self, scatter_df):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df,
            ax=ax,
            hue="cell_type",
            text_anno="cell_type",
            coding=False,
            show_legend=False,
        )
        assert isinstance(result, pd.DataFrame)

    def test_with_coding(self, scatter_df):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df,
            ax=ax,
            hue="cell_type",
            text_anno="cell_type",
            coding=True,
            show_legend=False,
        )
        assert "code" in result.columns

    def test_axis_format_empty(self, scatter_df):
        fig, ax = plt.subplots()
        categorical_scatter(
            scatter_df, ax=ax, hue="cell_type", axis_format="empty", show_legend=False
        )

    def test_no_hue(self, scatter_df):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df, ax=ax, color="steelblue", show_legend=False
        )
        assert isinstance(result, pd.DataFrame)

    def test_with_legend(self, scatter_df):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df,
            ax=ax,
            hue="cell_type",
            show_legend=True,
            legend_kws={"fontsize": 6},
        )
        assert isinstance(result, pd.DataFrame)

    def test_result_row_count(self, scatter_df):
        fig, ax = plt.subplots()
        result = categorical_scatter(
            scatter_df, ax=ax, hue="cell_type", show_legend=False
        )
        assert len(result) == len(scatter_df)


# ---------------------------------------------------------------------------
# continuous_scatter (utils)
# ---------------------------------------------------------------------------


class TestContinuousScatter:
    def test_runs_with_hue(self, scatter_df):
        fig, ax = plt.subplots()
        continuous_scatter(scatter_df, ax=ax, hue="n_genes")

    def test_runs_with_adata(self, adata):
        fig, ax = plt.subplots()
        continuous_scatter(adata, ax=ax, basis="umap", hue="n_genes")

    def test_no_hue(self, scatter_df):
        fig, ax = plt.subplots()
        continuous_scatter(scatter_df, ax=ax, color="steelblue")

    def test_custom_cmap(self, scatter_df):
        fig, ax = plt.subplots()
        continuous_scatter(scatter_df, ax=ax, hue="n_genes", cmap="plasma")

    def test_colorbar_false(self, scatter_df):
        fig, ax = plt.subplots()
        continuous_scatter(scatter_df, ax=ax, hue="n_genes", colorbar=False)


# ---------------------------------------------------------------------------
# plot_categorical (plotting)
# ---------------------------------------------------------------------------


class TestPlotCategorical:
    def test_runs(self, adata_with_colors):
        from adataviz.pl import plot_categorical

        plot_categorical(adata_with_colors, groupby="cell_type", show=False)

    def test_runs_with_explicit_palette(self, adata):
        from adataviz.pl import plot_categorical

        palette = {
            "TypeA": "#e41a1c",
            "TypeB": "#377eb8",
            "TypeC": "#4daf4a",
            "TypeD": "#984ea3",
            "TypeE": "#ff7f00",
        }
        # Pass palette via adata.uns so plot_categorical finds it
        adata.uns["cell_type_colors"] = list(palette.values())
        plot_categorical(adata, groupby="cell_type", show=False)

    def test_coding_false(self, adata_with_colors):
        from adataviz.pl import plot_categorical

        plot_categorical(
            adata_with_colors, groupby="cell_type", coding=False, show=False
        )

    def test_returns_axes(self, adata_with_colors):
        from adataviz.pl import plot_categorical
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        plot_categorical(adata_with_colors, ax=ax, groupby="cell_type", show=False)
        # Should not raise; ax is still a valid Axes object
        assert ax is not None


# ---------------------------------------------------------------------------
# Expression / heatmap plots (regressions for pandas>=3, category order, ...)
# ---------------------------------------------------------------------------

pch = pytest.importorskip("PyComplexHeatmap", reason="PyComplexHeatmap not installed")


@pytest.fixture
def adata_ordered(adata):
    """adata whose categorical order is deliberately *not* alphabetical."""
    cats = ["TypeE", "TypeC", "TypeA", "TypeD", "TypeB"]
    adata.obs["cell_type"] = adata.obs["cell_type"].cat.reorder_categories(cats)
    return adata


class TestExpressionPlots:
    genes = ["gene_0", "gene_1", "gene_2", "gene_3"]

    def test_complex_dotplot_runs_on_new_pandas(self, adata_ordered):
        # PyComplexHeatmap calls DataFrame.applymap, removed in pandas 3.
        from adataviz.pl import complex_dotplot

        cm = complex_dotplot(adata_ordered, self.genes, "cell_type", annotate=["condition"])
        assert cm.ax_heatmap is not None

    def test_gene_dotplot_runs(self, adata_ordered):
        from adataviz.pl import gene_dotplot

        cm = gene_dotplot(adata_ordered, self.genes, groupby="cell_type",
                          parent_col="condition", use_raw=False, title="t")
        assert cm.ax_heatmap is not None

    def test_complex_heatmap_keeps_category_order(self, adata_ordered):
        from adataviz.pl import complex_heatmap

        cm = complex_heatmap(adata_ordered, self.genes, "cell_type", row_cluster=False)
        assert list(cm.data2d.index) == list(adata_ordered.obs["cell_type"].cat.categories)

    def test_boxplot_keeps_category_order(self, adata_ordered):
        from adataviz.pl import boxplot

        ax = boxplot(adata_ordered, "gene_0", "cell_type", kind="box")
        labels = [t.get_text() for t in ax.get_xticklabels()]
        assert labels == list(adata_ordered.obs["cell_type"].cat.categories)

    def test_stacked_violin_both_orientations(self, adata_ordered):
        from adataviz.pl import stacked_violinplot

        axes = stacked_violinplot(adata_ordered, self.genes, "cell_type")
        assert len(axes) == len(self.genes)
        axes = stacked_violinplot(adata_ordered, self.genes, "cell_type", swap_axes=True)
        assert len(axes) == len(self.genes)

    def test_get_genes_mean_frac_matches_reference(self, adata_ordered):
        from adataviz.pl import get_genes_mean_frac

        out = get_genes_mean_frac(adata_ordered, self.genes, groupby="cell_type",
                                  use_raw=False, expression_cutoff=10)
        df = adata_ordered[:, self.genes].to_df()
        grp = adata_ordered.obs["cell_type"].astype(str)
        ref_mean = df.groupby(grp).mean()
        ref_frac = (df > 10).groupby(grp).mean()
        for _, row in out.iterrows():
            assert np.isclose(row["Mean"], ref_mean.loc[row["cell_type"], row["Gene"]])
            assert np.isclose(row["frac"], ref_frac.loc[row["cell_type"], row["Gene"]])
        # groups come back in categorical order, not alphabetically
        assert list(dict.fromkeys(out["cell_type"])) == list(
            adata_ordered.obs["cell_type"].cat.categories
        )

    def test_ridgeline_runs(self, adata_ordered):
        from adataviz.pl import ridgeline

        ax = ridgeline(adata_ordered, "gene_0", "cell_type")
        assert ax.figure is not None


class TestCompositionShow:
    def test_show_is_forwarded(self, adata_with_colors, monkeypatch):
        from adataviz.pl import area_plot, pie_plot, rose_plot, sankey_plot

        calls = []
        monkeypatch.setattr(plt, "show", lambda *a, **k: calls.append(1))
        rose_plot(adata_with_colors, "cell_type", show=True)
        pie_plot(adata_with_colors, "cell_type", split_by="condition", show=True)
        area_plot(adata_with_colors, "cell_type", "condition", show=True)
        sankey_plot(adata_with_colors, "cell_type", "condition", show=True)
        assert len(calls) == 4
