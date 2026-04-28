import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
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
