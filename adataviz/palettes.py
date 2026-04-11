from __future__ import annotations

import matplotlib.pylab as plt
import seaborn as sns
from matplotlib import cm, colors as mcolors
from matplotlib.colors import LinearSegmentedColormap, rgb2hex
import pandas as pd

# Colorblindness adjusted vega_10
# See https://github.com/scverse/scanpy/issues/387
vega_10 = list(map(mcolors.to_hex, cm.tab10.colors))
vega_10_scanpy = vega_10.copy()
vega_10_scanpy[2] = "#279e68"  # green
vega_10_scanpy[4] = "#aa40fc"  # purple
vega_10_scanpy[8] = "#b5bd61"  # kakhi

default_10 = vega_10_scanpy

# default matplotlib 2.0 palette
# see 'category20' on https://github.com/vega/vega/wiki/Scales#scale-range-literals
vega_20 = list(map(mcolors.to_hex, cm.tab20.colors))

# reorderd, some removed, some added
vega_20_scanpy = [
    # dark without grey:
    *vega_20[0:14:2],
    *vega_20[16::2],
    # light without grey:
    *vega_20[1:15:2],
    *vega_20[17::2],
    # manual additions:
    "#ad494a",
    "#8c6d31",
]
vega_20_scanpy[2] = vega_10_scanpy[2]
vega_20_scanpy[4] = vega_10_scanpy[4]
vega_20_scanpy[7] = vega_10_scanpy[8]  # kakhi shifted by missing grey

default_20 = vega_20_scanpy

# https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
# update 1
# orig reference https://research.wu.ac.at/en/publications/escaping-rgbland-selecting-colors-for-statistical-graphics-26
zeileis_28 = [
    "#023fa5",
    "#7d87b9",
    "#bec1d4",
    "#d6bcc0",
    "#bb7784",
    "#8e063b",
    "#4a6fe3",
    "#8595e1",
    "#b5bbe3",
    "#e6afb9",
    "#e07b91",
    "#d33f6a",
    "#11c638",
    "#8dd593",
    "#c6dec7",
    "#ead3c6",
    "#f0b98d",
    "#ef9708",
    "#0fcfc0",
    "#9cded6",
    "#d5eae7",
    "#f3e1eb",
    "#f6c4e1",
    "#f79cd4",
    # these last ones were added:
    "#7f7f7f",
    "#c7c7c7",
    "#1CE6FF",
    "#336600",
]

default_28 = zeileis_28

# from https://godsnotwheregodsnot.blogspot.com/2012/09/color-distribution-methodology.html
godsnot_102 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00",
    "#1CE6FF",
    "#FF34FF",
    "#FF4A46",
    "#008941",
    "#006FA6",
    "#A30059",
    "#FFDBE5",
    "#7A4900",
    "#0000A6",
    "#63FFAC",
    "#B79762",
    "#004D43",
    "#8FB0FF",
    "#997D87",
    "#5A0007",
    "#809693",
    "#6A3A4C",
    "#1B4400",
    "#4FC601",
    "#3B5DFF",
    "#4A3B53",
    "#FF2F80",
    "#61615A",
    "#BA0900",
    "#6B7900",
    "#00C2A0",
    "#FFAA92",
    "#FF90C9",
    "#B903AA",
    "#D16100",
    "#DDEFFF",
    "#000035",
    "#7B4F4B",
    "#A1C299",
    "#300018",
    "#0AA6D8",
    "#013349",
    "#00846F",
    "#372101",
    "#FFB500",
    "#C2FFED",
    "#A079BF",
    "#CC0744",
    "#C0B9B2",
    "#C2FF99",
    "#001E09",
    "#00489C",
    "#6F0062",
    "#0CBD66",
    "#EEC3FF",
    "#456D75",
    "#B77B68",
    "#7A87A1",
    "#788D66",
    "#885578",
    "#FAD09F",
    "#FF8A9A",
    "#D157A0",
    "#BEC459",
    "#456648",
    "#0086ED",
    "#886F4C",
    "#34362D",
    "#B4A8BD",
    "#00A6AA",
    "#452C2C",
    "#636375",
    "#A3C8C9",
    "#FF913F",
    "#938A81",
    "#575329",
    "#00FECF",
    "#B05B6F",
    "#8CD0FF",
    "#3B9700",
    "#04F757",
    "#C8A1A1",
    "#1E6E00",
    "#7900D7",
    "#A77500",
    "#6367A9",
    "#A05837",
    "#6B002C",
    "#772600",
    "#D790FF",
    "#9B9700",
    "#549E79",
    "#FFF69F",
    "#201625",
    "#72418F",
    "#BC23FF",
    "#99ADC0",
    "#3A2465",
    "#922329",
    "#5B4534",
    "#FDE8DC",
    "#404E55",
    "#0089A3",
    "#CB7E98",
    "#A4E804",
    "#324E72",
]

default_102 = godsnot_102


def create_spatial_colormap_viridis_style():
    """
    Creates a colormap similar to viridis but optimized for spatial transcriptomics.
    Background (0 expression) is dark blue/purple, increasing to bright yellow.
    """
    _colors = [
        "#0d0887",  # Dark purple (0 expression)
        "#2d1c8a",  # Deep purple
        "#46327e",  # Purple-blue
        "#5d4970",  # Mid purple
        "#756261",  # Purple-brown
        "#8e7d54",  # Brown-yellow
        "#a89947",  # Yellow-green
        "#c4b73b",  # Yellow
        "#e0d63b",  # Bright yellow
        "#fde724",  # Bright yellow (high expression)
    ]
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list("spatial_viridis", _colors, N=n_bins)
    return cmap


def create_spatial_colormap_plasma_style():
    """
    Creates a plasma-style colormap for spatial transcriptomics.
    Background is dark blue, increasing through purple to orange-red.
    """
    _colors = [
        "#0d0887",  # Dark blue (0 expression)
        "#350498",  # Deep purple
        "#5302a3",  # Purple
        "#7000a8",  # Magenta-purple
        "#8b0aa5",  # Magenta
        "#a31e9a",  # Pink-magenta
        "#b93289",  # Pink
        "#cc4678",  # Light pink
        "#db5c68",  # Salmon
        "#e97158",  # Orange-red
        "#f48849",  # Orange
        "#fba238",  # Yellow-orange
        "#febd2a",  # Yellow
        "#f0f921",  # Bright yellow (high expression)
    ]
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list("spatial_plasma", _colors, N=n_bins)
    return cmap


def create_spatial_colormap_fluorescence():
    """
    Creates a fluorescence microscopy-style colormap.
    Background is black, increasing through blue to cyan to white.
    """
    _colors = [
        "#000000",  # Black (0 expression)
        "#001122",  # Very dark blue
        "#002244",  # Dark blue
        "#003366",  # Deep blue
        "#004488",  # Blue
        "#0055aa",  # Bright blue
        "#0066cc",  # Cyan-blue
        "#0088ee",  # Light blue
        "#00aaff",  # Cyan
        "#44ccff",  # Light cyan
        "#88ddff",  # Very light cyan
        "#cceeff",  # Almost white cyan
        "#ffffff",  # White (high expression)
    ]
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list("spatial_fluor", _colors, N=n_bins)
    return cmap


def create_spatial_colormap_hot():
    """
    Creates a hot colormap with dark gray background.
    Background is dark gray, increasing through red to yellow to white.
    """
    _colors = [
        "#2a2a2a",  # Dark gray (0 expression)
        "#440000",  # Very dark red
        "#660000",  # Dark red
        "#880000",  # Red
        "#aa0000",  # Bright red
        "#cc2200",  # Red-orange
        "#ee4400",  # Orange
        "#ff6600",  # Bright orange
        "#ff8800",  # Yellow-orange
        "#ffaa00",  # Yellow
        "#ffcc00",  # Bright yellow
        "#ffee44",  # Light yellow
        "#ffffff",  # White (high expression)
    ]
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list("spatial_hot", _colors, N=n_bins)
    return cmap


def get_cmap(cmap):
    """
    Get a matplotlib colormap by name, compatible with different versions.

    Parameters
    ----------
    cmap : str
        Name of the colormap. Supports all standard matplotlib names
        ("viridis", "turbo", etc.) and custom registered colormaps
        ("parula", "exp1", "exp2", "meth1", "meth2").

    Returns
    -------
    matplotlib.colors.Colormap
        The requested colormap object.
    """
    try:
        return plt.colormaps.get(cmap)  # matplotlib >= 3.5.1?
    except AttributeError:
        return plt.get_cmap(cmap)  # matplotlib <=3.4.3?


def level_one_palette(name_list, order=None, palette="auto"):
    """
    Generate a color palette mapping for categorical names.

    Assigns a unique color to each unique value in ``name_list``
    using the specified palette.

    Parameters
    ----------
    name_list : array-like
        Series or list of category names (may contain duplicates/NaN).
    order : list, optional
        Desired order of categories. If None, categories are sorted
        alphabetically.
    palette : str or list, default "auto"
        Color palette to use. "auto" selects tab10/tab20/rainbow
        based on the number of unique categories. Can also be any
        seaborn/matplotlib palette name or a list of colors.

    Returns
    -------
    dict
        Mapping of category names to RGB color tuples.

    Raises
    ------
    ValueError
        If ``order`` doesn't match the unique values in ``name_list``.
    """
    name_set = set(name_list.dropna())
    if palette == "auto":
        if len(name_set) < 10:
            palette = "tab10"
        elif len(name_set) < 20:
            palette = "tab20"
        else:
            palette = "rainbow"

    if order is None:
        try:
            order = sorted(name_set)
        except TypeError:
            # name set contains multiple dtype (e.g., str and np.NaN)
            # do not sort
            order = list(name_set)
    else:
        if (set(order) != name_set) or (len(order) != len(name_set)):
            raise ValueError("Order is not equal to set(name_list).")

    n = len(order)
    _colors = sns.color_palette(palette, n)
    color_palette = {}
    for name, color in zip(order, _colors):
        color_palette[name] = color
    return color_palette


def add_color_scheme(
    adata,
    column: str,
    palette_key: str | None = None,
    palette: str | list | dict | None = None,
    na_color: str = "#CCCCCC",
):
    """
    Add a color palette to AnnData.uns based on an existing categorical column.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    column : str
        Name of the column in adata.obs to create colors for
    palette_key : str, optional
        Key name for the palette in adata.uns. If None, uses f"{column}_colors"
    palette : str, list, or dict, optional
        Color palette to use. Can be:
        - matplotlib/seaborn palette name (e.g., 'tab10', 'Set1', 'viridis')
        - list of hex colors
        - dict mapping category values to hex colors
        - None (uses default palettes sized to the number of categories)
    na_color : str, default '#CCCCCC'
        Color to use for NaN/missing values and for categories not found in a
        dict palette

    Returns
    -------
    dict or str
        For categorical data: dict mapping category → hex color, including
        NaN → na_color if missing values are present.
        For continuous data: the colormap name or palette passed in.

    Examples
    --------
    >>> add_color_scheme(adata, 'cell_type')
    >>> add_color_scheme(adata, 'cluster', palette='Set1')
    >>> add_color_scheme(adata, 'condition', palette={'A': '#FF0000', 'B': '#0000FF'})
    """
    if column not in adata.obs.columns:
        raise ValueError(f"Column '{column}' not found in adata.obs")

    if palette_key is None:
        palette_key = f"{column}_colors"

    if adata.uns is None:
        adata.uns = {}

    data = adata.obs[column].copy()
    has_na = data.isna().any()

    if (
        isinstance(data.dtype, pd.CategoricalDtype)
        or pd.api.types.is_object_dtype(data)
        or isinstance(data.dtype, pd.StringDtype)
    ):
        if isinstance(data.dtype, pd.CategoricalDtype):
            unique_vals = data.cat.categories.tolist()
        else:
            unique_vals = data.dropna().unique()

        n_categories = len(unique_vals)

        if isinstance(palette, dict):
            color_map = {cat: palette.get(cat, na_color) for cat in unique_vals}
            if has_na:
                color_map[pd.NA] = na_color
            adata.uns[palette_key] = [color_map[cat] for cat in unique_vals]
            return color_map

        if palette is None:
            if n_categories <= 10:
                palette = default_10
            elif n_categories <= 20:
                palette = default_20
            elif n_categories <= 28:
                palette = default_28
            else:
                palette = default_102
            _colors = sns.color_palette(palette, n_categories)
        elif isinstance(palette, str):
            _colors = sns.color_palette(palette, n_categories)
        elif isinstance(palette, list):
            if len(palette) < n_categories:
                _colors = (palette * ((n_categories // len(palette)) + 1))[
                    :n_categories
                ]
            else:
                _colors = palette[:n_categories]
        else:
            raise ValueError("palette must be a string, list, or dict of colors")

        hex_colors = [rgb2hex(c) if not isinstance(c, str) else c for c in _colors]
        adata.uns[palette_key] = hex_colors

        result = dict(zip(unique_vals, hex_colors))
        if has_na:
            result[pd.NA] = na_color
        return result

    else:
        if palette is None:
            palette = "viridis"
        adata.uns[palette_key] = palette
        return palette


def add_colors(adata, cat_col, palette):
    """
    Assign colors to AnnData categorical column in ``uns``.

    Maps each category in ``adata.obs[cat_col]`` to a color from
    ``palette`` and stores the list in ``adata.uns[f"{cat_col}_colors"]``.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with a categorical column in ``.obs``.
    cat_col : str
        Name of the categorical column in ``adata.obs``.
    palette : dict or pd.DataFrame
        Color mapping. If dict, maps category names to hex colors.
        If DataFrame, must have a "Hex" column with category names
        as index.
    """
    _colors = []
    for _cat in adata.obs[cat_col].cat.categories:
        try:
            if isinstance(palette, dict):
                color = palette[_cat]
            else:
                color = palette.loc[_cat, "Hex"]
        except KeyError:
            print(_cat)
            color = "#808080"
        _colors.append(color)

    adata.uns[f"{cat_col}_colors"] = _colors
