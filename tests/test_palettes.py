import re

import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import LinearSegmentedColormap, is_color_like

from adataviz import palettes


# ---------------------------------------------------------------------------
# Palette constants
# ---------------------------------------------------------------------------

HEX_RE = re.compile(r"^#[0-9a-fA-F]{6}$")


def is_hex(s):
    return bool(HEX_RE.match(s))


def test_default_10_length():
    assert len(palettes.default_10) == 10


def test_default_20_length():
    assert len(palettes.default_20) == 20


def test_default_28_length():
    assert len(palettes.default_28) == 28


def test_default_102_length():
    assert len(palettes.default_102) == 102


def test_palette_constants_are_valid_colors():
    for palette in (
        palettes.default_10,
        palettes.default_20,
        palettes.default_28,
        palettes.default_102,
    ):
        for c in palette:
            assert is_color_like(c), f"Not a valid color: {c}"


def test_vega_10_scanpy_overrides():
    # Colorblindness-adjusted overrides should be present
    assert palettes.vega_10_scanpy[2] == "#279e68"
    assert palettes.vega_10_scanpy[4] == "#aa40fc"
    assert palettes.vega_10_scanpy[8] == "#b5bd61"


# ---------------------------------------------------------------------------
# Spatial colormap factories
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory",
    [
        palettes.create_spatial_colormap_viridis_style,
        palettes.create_spatial_colormap_plasma_style,
        palettes.create_spatial_colormap_fluorescence,
        palettes.create_spatial_colormap_hot,
    ],
)
def test_spatial_colormaps_return_lsc(factory):
    cmap = factory()
    assert isinstance(cmap, LinearSegmentedColormap)


def test_spatial_colormaps_have_256_bins():
    for factory in (
        palettes.create_spatial_colormap_viridis_style,
        palettes.create_spatial_colormap_plasma_style,
        palettes.create_spatial_colormap_fluorescence,
        palettes.create_spatial_colormap_hot,
    ):
        cmap = factory()
        assert cmap.N == 256


# ---------------------------------------------------------------------------
# get_cmap
# ---------------------------------------------------------------------------


def test_get_cmap_returns_colormap():
    cmap = palettes.get_cmap("viridis")
    assert cmap is not None
    # calling it with a scalar should return an RGBA tuple
    result = cmap(0.5)
    assert len(result) == 4


def test_get_cmap_unknown_returns_none():
    # plt.colormaps.get returns None for unknown names rather than raising
    result = palettes.get_cmap("this_cmap_does_not_exist_xyz")
    assert result is None


# ---------------------------------------------------------------------------
# level_one_palette
# ---------------------------------------------------------------------------


def test_level_one_palette_returns_dict(scatter_df):
    result = palettes.level_one_palette(scatter_df["cell_type"])
    assert isinstance(result, dict)
    assert set(result.keys()) == set(scatter_df["cell_type"].dropna().unique())


def test_level_one_palette_auto_small(scatter_df):
    # 5 categories → should pick tab10 (len < 10)
    result = palettes.level_one_palette(scatter_df["cell_type"])
    assert len(result) == 5


def test_level_one_palette_explicit_palette(scatter_df):
    result = palettes.level_one_palette(scatter_df["cell_type"], palette="Set2")
    assert len(result) == 5
    for v in result.values():
        assert is_color_like(v)


def test_level_one_palette_explicit_order(scatter_df):
    order = ["TypeA", "TypeB", "TypeC", "TypeD", "TypeE"]
    result = palettes.level_one_palette(scatter_df["cell_type"], order=order)
    assert list(result.keys()) == order


def test_level_one_palette_invalid_order_raises(scatter_df):
    with pytest.raises(ValueError, match="Order is not equal"):
        palettes.level_one_palette(scatter_df["cell_type"], order=["TypeA", "TypeB"])


# ---------------------------------------------------------------------------
# add_color_scheme — categorical column
# ---------------------------------------------------------------------------


def test_add_color_scheme_default_palette(adata):
    result = palettes.add_color_scheme(adata, "cell_type")
    assert isinstance(result, dict)
    assert set(result.keys()) == set(adata.obs["cell_type"].cat.categories)
    assert "cell_type_colors" in adata.uns
    assert len(adata.uns["cell_type_colors"]) == 5


def test_add_color_scheme_str_palette(adata):
    result = palettes.add_color_scheme(adata, "cell_type", palette="Set1")
    assert len(result) == 5
    for v in result.values():
        assert is_color_like(v)


def test_add_color_scheme_list_palette(adata):
    palette = ["#ff0000", "#00ff00", "#0000ff", "#ffff00", "#ff00ff"]
    result = palettes.add_color_scheme(adata, "cell_type", palette=palette)
    assert list(result.values()) == palette


def test_add_color_scheme_list_too_short(adata):
    # List shorter than n_categories → should tile and truncate
    palette = ["#ff0000", "#00ff00"]
    result = palettes.add_color_scheme(adata, "cell_type", palette=palette)
    assert len(result) == 5


def test_add_color_scheme_dict_palette(adata):
    palette = {
        "TypeA": "#ff0000",
        "TypeB": "#00ff00",
        "TypeC": "#0000ff",
        "TypeD": "#ffff00",
        "TypeE": "#ff00ff",
    }
    result = palettes.add_color_scheme(adata, "cell_type", palette=palette)
    assert result["TypeA"] == "#ff0000"
    assert result["TypeE"] == "#ff00ff"


def test_add_color_scheme_dict_missing_category_gets_na_color(adata):
    # Only map 3 out of 5 categories
    palette = {"TypeA": "#ff0000", "TypeB": "#00ff00", "TypeC": "#0000ff"}
    na_color = "#AABBCC"
    result = palettes.add_color_scheme(
        adata, "cell_type", palette=palette, na_color=na_color
    )
    assert result["TypeD"] == na_color
    assert result["TypeE"] == na_color


def test_add_color_scheme_na_in_data(adata):
    # Insert NaN into a non-categorical object column
    adata.obs["cell_type_obj"] = adata.obs["cell_type"].astype(str)
    adata.obs.loc["cell_0", "cell_type_obj"] = np.nan
    na_color = "#CCCCCC"
    result = palettes.add_color_scheme(adata, "cell_type_obj", na_color=na_color)
    # The NaN key should be present in the returned dict
    nan_keys = [k for k in result.keys() if k is pd.NA]
    assert len(nan_keys) == 1
    assert result[nan_keys[0]] == na_color


def test_add_color_scheme_dict_na_in_data(adata):
    # Dict palette + NaN values → both unmatched and NaN get na_color
    adata.obs["cell_type_obj"] = adata.obs["cell_type"].astype(str)
    adata.obs.loc["cell_0", "cell_type_obj"] = np.nan
    palette = {"TypeA": "#ff0000"}
    na_color = "#DDDDDD"
    result = palettes.add_color_scheme(
        adata, "cell_type_obj", palette=palette, na_color=na_color
    )
    nan_keys = [k for k in result.keys() if k is pd.NA]
    assert len(nan_keys) == 1
    assert result[nan_keys[0]] == na_color


def test_add_color_scheme_stores_in_uns(adata):
    custom_key = "my_colors"
    palettes.add_color_scheme(adata, "cell_type", palette_key=custom_key)
    assert custom_key in adata.uns


def test_add_color_scheme_invalid_column(adata):
    with pytest.raises(ValueError, match="not found in adata.obs"):
        palettes.add_color_scheme(adata, "nonexistent_column")


def test_add_color_scheme_invalid_palette_type(adata):
    with pytest.raises(ValueError, match="palette must be"):
        palettes.add_color_scheme(adata, "cell_type", palette=12345)


# ---------------------------------------------------------------------------
# add_color_scheme — continuous column
# ---------------------------------------------------------------------------


def test_add_color_scheme_continuous_default(adata):
    result = palettes.add_color_scheme(adata, "n_genes")
    assert result == "viridis"
    assert adata.uns["n_genes_colors"] == "viridis"


def test_add_color_scheme_continuous_custom(adata):
    result = palettes.add_color_scheme(adata, "n_genes", palette="plasma")
    assert result == "plasma"


# ---------------------------------------------------------------------------
# add_colors
# ---------------------------------------------------------------------------


def test_add_colors_dict(adata):
    palette = {
        cat: f"#{i:02x}{i:02x}{i:02x}"
        for i, cat in enumerate(adata.obs["cell_type"].cat.categories)
    }
    palettes.add_colors(adata, "cell_type", palette)
    assert "cell_type_colors" in adata.uns
    assert len(adata.uns["cell_type_colors"]) == 5


def test_add_colors_missing_key_gets_gray(adata):
    # Only partial dict → missing categories fall back to #808080
    palette = {"TypeA": "#ff0000"}
    palettes.add_colors(adata, "cell_type", palette)
    stored = adata.uns["cell_type_colors"]
    # TypeB, TypeC, TypeD, TypeE are missing → gray
    assert stored[1] == "#808080"
