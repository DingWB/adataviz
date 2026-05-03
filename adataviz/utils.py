import pandas as pd
import os
import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from matplotlib.text import Text
import copy
import warnings

# Note: do NOT call ``warnings.filterwarnings("ignore")`` at module import.
# Library code must not silence warnings globally for downstream users.
# Use ``with warnings.catch_warnings():`` locally where needed.

mm2inch = 1 / 25.4


def df2stdout(df):
    """
    Write a DataFrame to stdout in tab-separated format.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write. Each row is printed as a tab-separated line.
        NaN values are replaced with empty strings.
    """
    sys.stdout.write("\t".join([str(i) for i in df.columns.tolist()]) + "\n")
    for i, row in df.iterrows():
        try:
            sys.stdout.write(
                "\t".join([str(i) for i in row.fillna("").tolist()]) + "\n"
            )
        except (TypeError, UnicodeEncodeError, OSError, ValueError):
            sys.stdout.close()


def serialize(x):
    """
    Serialize an object to stdout.

    If ``x`` is a DataFrame, writes it as tab-separated text via
    :func:`df2stdout`. Otherwise prints ``x`` directly. Does nothing
    if ``x`` is None.

    Parameters
    ----------
    x : pd.DataFrame, any, or None
        Object to serialize.
    """
    if isinstance(x, pd.DataFrame):
        df2stdout(x)
    elif x is not None:
        print(x)
    else:
        pass


def prepare_color_palette(color_dict=None, outpath="palette.xlsx"):
    """
    Generate an Excel file with colored color palette sheets.

    Each key in ``color_dict`` becomes a sheet, with category names as
    the index and hex color codes in the "Hex" column. Cells are
    formatted with their corresponding background color for visual
    reference.

    Parameters
    ----------
    color_dict : dict of dict
        Nested dictionary where outer keys are sheet names (e.g.
        "Class", "Subclass") and inner dicts map category names to
        hex color codes, e.g.
        ``{"Class": {"Excitatory": "#FF0000", "Inhibitory": "#0000FF"}}``.
    outpath : str, default "palette.xlsx"
        Output path for the Excel file. Supports ``~`` expansion.
    """
    outpath = os.path.expanduser(outpath)
    writer = pd.ExcelWriter(outpath)
    for key in color_dict:
        data = pd.DataFrame.from_dict(color_dict[key], orient="index", columns=["Hex"])
        # data.style.background_gradient(cmap='gray_r')
        # data.style.applymap(lambda x:'color:'+x if x.startswith('#') else 'color: white')
        data.to_excel(writer, sheet_name=key, index=True)
        workbook = writer.book
        worksheet = writer.sheets[key]
        colors = data.Hex.tolist()
        for i in range(data.shape[0]):
            color = colors[i]
            f = workbook.add_format(
                {"bold": True, "font_color": "black", "bg_color": color}
            )
            worksheet.write(i + 1, 1, color, f)
        width = 20
        cell_fmt = workbook.add_format(
            {
                "bold": False,
                "font_color": "black",
                # 'bg_color':'green',
                "align": "center",
                "valign": "vcenter",
            }
        )
        # styled = data.style.applymap(lambda val: 'color: %s' % 'red' if val < 0 else 'black').highlight_max()
        worksheet.set_column(0, 1, width, cell_fmt)
    # worksheet.conditional_format(f'A:{last_col}', {'type': 'no_blanks', 'format': cell_fmt})
    writer.close()


def mpl_style():
    """
    Set publication-quality matplotlib defaults.

    Configures matplotlib to use:

    - Font type 42 for PDF/PS (editable text in Illustrator)
    - Arial sans-serif font
    - 80 DPI for display, 300 DPI for saved figures
    """
    import matplotlib as mpl

    mpl.style.use("default")
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = "Arial"
    mpl.rcParams["figure.dpi"] = 80
    mpl.rcParams["savefig.dpi"] = 300


def parse_json(url):
    """
    Parse Allen Brain Atlas JSON ontology into a DataFrame.

    Parameters
    ----------
    url : str
        URL to the Allen Brain Atlas ontology JSON endpoint,
        e.g. ``"https://atlas.brain-map.org/atlasviewer/ontologies/11.json"``.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: acronym, Hex, id, name,
        parent_structure_id, safe_name, structure_id_path.
    """
    import json
    import requests

    data = json.loads(requests.get(url).content.decode())
    R = []
    for record in data["msg"]:
        R.append(
            [
                record["acronym"],
                record["color_hex_triplet"],
                record["id"],
                record["name"],
                record["parent_structure_id"],
                record["safe_name"],
                record["structure_id_path"],
            ]
        )
    df = pd.DataFrame(
        R,
        columns=[
            "acronym",
            "Hex",
            "id",
            "name",
            "parent_structure_id",
            "safe_name",
            "structure_id_path",
        ],
    )
    return df


def get_brain_region_structure():
    """
    https://atlas.brain-map.org/
    BICAN: https://atlas.brain-map.org/atlasviewer/ontologies/11.json
    HBA: https://atlas.brain-map.org/atlasviewer/ontologies/7.json
    Returns
    -------

    """
    # pip install XlsxWriter
    writer = pd.ExcelWriter("AllenBrainRegionStructure.xlsx")
    for url, name in zip(
        [
            "https://atlas.brain-map.org/atlasviewer/ontologies/11.json",
            "https://atlas.brain-map.org/atlasviewer/ontologies/7.json",
        ],
        ["BICAN_Brodmann", "HBA_guide"],
    ):
        print(name)
        df = parse_json(url)
        df["Hex"] = "#" + df.Hex.map(str)
        df["id"] = df["id"].fillna("None").map(str)
        id2acronym = df.set_index("id").acronym.to_dict()
        df.insert(
            0,
            "Parent",
            df.parent_structure_id.fillna(0).map(int).map(str).map(id2acronym),
        )
        df.parent_structure_id = df.parent_structure_id.map(str)
        df.structure_id_path = df.structure_id_path.apply(lambda x: x[1:-1])
        df["structure_path"] = df.structure_id_path.apply(
            lambda x: "//".join([id2acronym[p] for p in x.split("/") if p != ""])
        )
        # data.style.background_gradient(cmap='gray_r')
        # data.style.applymap(lambda x:'color:'+x if x.startswith('#') else 'color: white')
        df = df.loc[
            :,
            [
                "Parent",
                "acronym",
                "Hex",
                "name",
                "safe_name",
                "id",
                "parent_structure_id",
                "structure_path",
                "structure_id_path",
            ],
        ]
        df.to_excel(writer, sheet_name=name, index=False)
        workbook = writer.book
        worksheet = writer.sheets[name]
        colors = df.Hex.tolist()
        col_idx = df.columns.tolist().index("Hex")
        for i in range(df.shape[0]):
            color = colors[i]
            f = workbook.add_format(
                {"bold": True, "font_color": "black", "bg_color": color}
            )
            worksheet.write(
                i + 1, col_idx, color, f
            )  # worksheet.write(row, col, *args), row and col are 0-based
        # width=20
        # cell_fmt = workbook.add_format(
        #                 {'bold': False,'font_color': 'black',
        #                 # 'bg_color':'green',
        #                 'align': 'center', 'valign': 'vcenter'})
        # # styled = data.style.applymap(lambda val: 'color: %s' % 'red' if val < 0 else 'black').highlight_max()
        # worksheet.set_column(0,col_idx,width,cell_fmt) #first_col,last_col
        # worksheet.conditional_format(f'A:{last_col}', {'type': 'no_blanks', 'format': cell_fmt})
    writer.close()
    df_bican = pd.read_excel("Jon.xlsx")
    df_bican.Acronym = df_bican.Acronym.apply(lambda x: x.split("(")[0].strip())
    regions = df_bican.Acronym.tolist()
    df = pd.read_excel("AllenBrainRegionStructure.xlsx", sheet_name="BICAN_Brodmann")
    # df['Keep']=df.structure_path.apply(lambda x:[r for r in x.split('//') if r in regions])
    # df=df.loc[df.Keep.apply(len) > 0]
    # df.drop('Keep',axis=1,inplace=True)
    df = df.loc[df.acronym.isin(regions)]
    df.parent_structure_id = df.parent_structure_id.map(int)

    writer = pd.ExcelWriter("BICAN_regions.xlsx")
    df.to_excel(writer, sheet_name="BICAN", index=False)
    workbook = writer.book
    worksheet = writer.sheets["BICAN"]
    colors = df.Hex.tolist()
    col_idx = df.columns.tolist().index("Hex")
    for i in range(df.shape[0]):
        color = colors[i]
        f = workbook.add_format(
            {"bold": True, "font_color": "black", "bg_color": color}
        )
        worksheet.write(i + 1, col_idx, color, f)
    writer.close()


def read_google_sheet(url=None, **kwargs):
    """
    Read a Google Sheets tab as a pandas DataFrame.

    Exports the sheet as TSV and reads it via ``pd.read_csv``.

    Parameters
    ----------
    url : str
        Full Google Sheets URL including the ``gid`` parameter,
        e.g. ``"https://docs.google.com/spreadsheets/d/<ID>/edit?gid=<GID>#gid=<GID>"``.
    **kwargs
        Additional keyword arguments passed to ``pd.read_csv``.

    Returns
    -------
    pd.DataFrame
        Contents of the specified sheet tab.
    """
    assert url is not None
    # url="https://docs.google.com/spreadsheets/d/12H3p2F_qrcQ3ymF614VRzU_6vcVTsRca0uOIBQJXVaU/edit?gid=1969763406#gid=1969763406"
    Id = url.split("/d/")[1].split("/")[0]
    gid = url.split("?gid=")[1].split("#gid=")[0]
    df = pd.read_csv(
        f"https://docs.google.com/spreadsheets/d/{Id}/export?format=tsv&id={Id}&gid={gid}",
        sep="\t",
        **kwargs,
    )
    return df


def despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False):
    """
    Remove the top and right spines from plot(s).

    Parameters
    ----------
    fig : matplotlib figure, optional
            Figure to despine all axes of, defaults to the current figure.
    ax : matplotlib axes, optional
            Specific axes object to despine. Ignored if fig is provided.
    top, right, left, bottom : boolean, optional
            If True, remove that spine.

    Returns
    -------
    None

    """
    if fig is None and ax is None:
        axes = plt.gcf().axes
    elif fig is not None:
        axes = fig.axes
    elif ax is not None:
        axes = [ax]

    for ax_i in axes:
        for side in ["top", "right", "left", "bottom"]:
            is_visible = not locals()[side]
            ax_i.spines[side].set_visible(is_visible)
        if left and not right:  # remove left, keep right
            maj_on = any(t.tick1line.get_visible() for t in ax_i.yaxis.majorTicks)
            min_on = any(t.tick1line.get_visible() for t in ax_i.yaxis.minorTicks)
            ax_i.yaxis.set_ticks_position("right")
            for t in ax_i.yaxis.majorTicks:
                t.tick2line.set_visible(maj_on)
            for t in ax_i.yaxis.minorTicks:
                t.tick2line.set_visible(min_on)

        if bottom and not top:
            maj_on = any(t.tick1line.get_visible() for t in ax_i.xaxis.majorTicks)
            min_on = any(t.tick1line.get_visible() for t in ax_i.xaxis.minorTicks)
            ax_i.xaxis.set_ticks_position("top")
            for t in ax_i.xaxis.majorTicks:
                t.tick2line.set_visible(maj_on)
            for t in ax_i.xaxis.minorTicks:
                t.tick2line.set_visible(min_on)


def _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=5):
    """
    Add small arrow axis labels to the bottom-left corner of a scatter plot.

    Removes default axis ticks/spines and draws two small arrows with labels
    indicating the x and y coordinate names (e.g. "UMAP 1", "UMAP 2").

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to annotate.
    x : str
        Label for the horizontal arrow (typically the basis name like "umap_0").
        Underscores are replaced with spaces and text is uppercased.
    y : str
        Label for the vertical arrow.
    arrow_kws : dict, optional
        Keyword arguments passed to ``ax.arrow()``, overriding defaults
        (width=0.003, linewidth=0, color="black").
    fontsize : float, default 5
        Font size for the axis labels.
    """
    # clean ax axises
    ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    despine(ax=ax, left=True, bottom=True)

    _arrow_kws = {"width": 0.003, "linewidth": 0, "color": "black"}
    if arrow_kws is not None:
        _arrow_kws.update(arrow_kws)

    ax.arrow(0.06, 0.06, 0, 0.06, **_arrow_kws, transform=ax.transAxes)
    ax.arrow(0.06, 0.06, 0.06, 0, **_arrow_kws, transform=ax.transAxes)
    ax.text(
        0.06,
        0.03,
        x.upper().replace("_", " "),
        fontdict={
            "fontsize": fontsize,
            "horizontalalignment": "left",
            "verticalalignment": "center",
        },
        transform=ax.transAxes,
    )
    ax.text(
        0.03,
        0.06,
        y.upper().replace("_", " "),
        fontdict={
            "fontsize": fontsize,
            "rotation": 90,
            "rotation_mode": "anchor",
            "horizontalalignment": "left",
            "verticalalignment": "center",
        },
        transform=ax.transAxes,
    )
    return


def zoom_min_max(vmin, vmax, scale):
    """
    Symmetrically expand a value range by a scale factor.

    Parameters
    ----------
    vmin : float
        Minimum value of the range.
    vmax : float
        Maximum value of the range.
    scale : float
        Scale factor. 1.0 = no change, 1.2 = 20% wider, etc.

    Returns
    -------
    tuple of (float, float)
        The expanded (vmin, vmax) range.
    """
    width = vmax - vmin
    width_zoomed = width * scale
    delta_value = (width_zoomed - width) / 2
    return vmin - delta_value, vmax + delta_value


def zoom_ax(ax, zoom_scale, on="both"):
    """
    Zoom (expand) the axis limits of a matplotlib Axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to zoom.
    zoom_scale : float
        Scale factor for the axis limits. Values > 1 zoom out
        (show more space), values < 1 zoom in.
    on : str, default "both"
        Which axes to zoom. Options: "both", "x", "y".
    """
    on = on.lower()
    xlim = ax.get_xlim()
    xlim_zoomed = zoom_min_max(vmin=xlim[0], vmax=xlim[1], scale=zoom_scale)

    ylim = ax.get_ylim()
    ylim_zoomed = zoom_min_max(vmin=ylim[0], vmax=ylim[1], scale=zoom_scale)

    if (on == "both") or ("x" in on):
        ax.set_xlim(xlim_zoomed)
    if (on == "both") or ("y" in on):
        ax.set_ylim(ylim_zoomed)


def _extract_coords(data, coord_base, x, y):
    """
    Extract 2D coordinates from various data formats.

    Supports AnnData (via ``obsm["X_{coord_base}"]``), xarray Dataset,
    and plain DataFrames.

    Parameters
    ----------
    data : anndata.AnnData, xr.Dataset, or pd.DataFrame
        Input data containing coordinates.
    coord_base : str
        Name of the coordinate embedding, e.g. "umap", "tsne".
        For AnnData, looks up ``obsm["X_{coord_base}"]``.
    x : str or None
        Column name for x coordinate. If None, defaults to
        ``"{coord_base}_0"``.
    y : str or None
        Column name for y coordinate. If None, defaults to
        ``"{coord_base}_1"``.

    Returns
    -------
    tuple of (pd.DataFrame, str, str)
        - DataFrame with columns "x" and "y"
        - x column name used
        - y column name used
    """
    import xarray as xr
    import anndata

    if (x is not None) and (y is not None):
        pass
    else:
        x = f"{coord_base}_0"
        y = f"{coord_base}_1"

    if isinstance(data, anndata.AnnData):
        adata = data
        _data = pd.DataFrame(
            {
                "x": adata.obsm[f"X_{coord_base}"][:, 0],
                "y": adata.obsm[f"X_{coord_base}"][:, 1],
            },
            index=adata.obs_names,
        )
    elif isinstance(data, xr.Dataset):
        ds = data
        if coord_base not in ds.dims:
            raise KeyError(f"xr.Dataset do not contain {coord_base} dim")
        data_var = {i for i in ds.data_vars.keys() if i.startswith(coord_base)}.pop()
        _data = pd.DataFrame(
            {
                "x": ds[data_var].sel({coord_base: f"{coord_base}_0"}).to_pandas(),
                "y": ds[data_var].sel({coord_base: f"{coord_base}_1"}).to_pandas(),
            }
        )
    else:
        if (x not in data.columns) or (y not in data.columns):
            raise KeyError(f"{x} or {y} not found in columns.")
        _data = pd.DataFrame({"x": data[x], "y": data[y]})
    return _data, x, y


def _density_based_sample(
    data: pd.DataFrame, coords: list, portion=None, size=None, seed=None
):
    """
    Downsample data based on local density to reduce overplotting.

    Points in dense regions have lower probability of being selected,
    while isolated points are more likely to be kept. Uses Local
    Outlier Factor (LOF) to estimate density.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing at least the coordinate columns.
    coords : list of str
        Column names to use as coordinates for density estimation.
    portion : float, optional
        Fraction of data to keep (0 to 1). Mutually exclusive with
        ``size``.
    size : int, optional
        Exact number of points to keep. Mutually exclusive with
        ``portion``.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Downsampled DataFrame with density-weighted selection.

    Raises
    ------
    ValueError
        If neither ``portion`` nor ``size`` is provided.
    """
    from sklearn.neighbors import LocalOutlierFactor

    clf = LocalOutlierFactor(
        n_neighbors=20,
        algorithm="auto",
        leaf_size=30,
        metric="minkowski",
        p=2,
        metric_params=None,
        contamination=0.1,
    )

    # coords should already exist in data, get them by column names list
    data_coords = data[coords]
    clf.fit(data_coords)
    # original score is negative, the larger the denser
    density_score = clf.negative_outlier_factor_
    delta = density_score.max() - density_score.min()
    # density score to probability: the denser the less probability to be picked up
    probability_score = 1 - (density_score - density_score.min()) / delta
    probability_score = np.sqrt(probability_score)
    probability_score = probability_score / probability_score.sum()

    if size is not None:
        pass
    elif portion is not None:
        size = int(data_coords.index.size * portion)
    else:
        raise ValueError("Either portion or size should be provided.")
    if seed is not None:
        np.random.seed(seed)
    selected_cell_index = np.random.choice(
        data_coords.index, size=size, replace=False, p=probability_score
    )  # choice data based on density weights

    # return the down sampled data
    return data.reindex(selected_cell_index)


def _auto_size(ax, n_dots):
    """
    Automatically determine scatter dot size based on axes dimensions and point count.

    Scales dot size inversely with the number of points, adjusted for
    the physical size of the axes. Larger figures or fewer points
    result in bigger dots.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to measure physical size from.
    n_dots : int
        Number of dots to be plotted.

    Returns
    -------
    float
        Recommended dot size (``s`` parameter for scatter plots).
    """
    bbox = ax.get_window_extent().transformed(
        ax.get_figure().dpi_scale_trans.inverted()
    )
    scale = bbox.width * bbox.height / 14.6  # 14.6 is a 5*5 fig I used to estimate
    n = n_dots / scale  # larger figure means data look sparser
    if n < 500:
        s = 14 - n / 100
    elif n < 1500:
        s = 7
    elif n < 3000:
        s = 5
    elif n < 8000:
        s = 3
    elif n < 15000:
        s = 2
    elif n < 30000:
        s = 1.5
    elif n < 50000:
        s = 1
    elif n < 80000:
        s = 0.8
    elif n < 150000:
        s = 0.6
    elif n < 300000:
        s = 0.5
    elif n < 500000:
        s = 0.4
    elif n < 800000:
        s = 0.3
    elif n < 1000000:
        s = 0.2
    elif n < 2000000:
        s = 0.1
    elif n < 3000000:
        s = 0.07
    elif n < 4000000:
        s = 0.05
    elif n < 5000000:
        s = 0.03
    else:
        s = 0.02
    return s


def _take_data_series(data, k):
    """
    Extract a named series from various data formats.

    Parameters
    ----------
    data : anndata.AnnData, xr.Dataset, xr.DataArray, or pd.DataFrame
        Input data object.
    k : str
        Key/column name to extract. For AnnData, looks up ``obs[k]``.

    Returns
    -------
    pd.Series
        Copy of the requested data series.
    """
    import xarray as xr
    import anndata

    if isinstance(data, (xr.Dataset, xr.DataArray)):
        _value = data[k].to_pandas()
    elif isinstance(data, anndata.AnnData):
        _value = data.obs[k].copy()
    else:
        _value = data[k].copy()
    return _value


def level_one_palette(name_list, order=None, palette="auto"):
    """
    Generate a color palette for a list of category names.

    .. deprecated::
        Use :func:`adataviz.palettes.level_one_palette` instead.

    Parameters
    ----------
    name_list : array-like
        List or array of category names.
    order : list, optional
        Desired order of categories. If None, sorted alphabetically.
    palette : str or list, default "auto"
        Palette name or list of colors. "auto" selects from
        scanpy-style palettes based on the number of categories.

    Returns
    -------
    dict
        Mapping of category names to hex color strings.
    """
    from . import palettes as _palettes

    warnings.warn(
        "adataviz.utils.level_one_palette is deprecated; use adataviz.palettes.level_one_palette instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _palettes.level_one_palette(name_list, order=order, palette=palette)


def _calculate_luminance(color):
    """
    Calculate the relative luminance of a color according to W3C standards

    Parameters
    ----------
    color : matplotlib color or sequence of matplotlib colors
            Hex code, rgb-tuple, or html color name.
    Returns
    -------
    luminance : float(s) between 0 and 1

    """
    rgb = matplotlib.colors.colorConverter.to_rgba_array(color)[:, :3]
    rgb = np.where(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    lum = rgb.dot([0.2126, 0.7152, 0.0722])
    try:
        return lum.item()
    except ValueError:
        return lum


def _text_anno_scatter(
    data: pd.DataFrame,
    ax,
    x: str,
    y: str,
    dodge_text=False,
    anno_col="text_anno",
    text_kws=None,
    text_transform=None,
    dodge_kws=None,
    luminance=0.48,
):
    """Add text annotation to a scatter plot."""
    # prepare kws
    text_kws = {} if text_kws is None else text_kws
    text_kws.setdefault("fontsize", 5)
    text_kws.setdefault("fontweight", "black")
    text_kws.setdefault("ha", "center")  # horizontalalignment
    text_kws.setdefault("va", "center")  # verticalalignment
    text_kws.setdefault("color", "black")  # c
    bbox = dict(
        boxstyle="round",
        edgecolor=(0.5, 0.5, 0.5, 0.2),
        fill=False,
        facecolor=(0.8, 0.8, 0.8, 0.2),
        alpha=1,
        linewidth=0.5,
    )
    text_kws.setdefault("bbox", bbox)
    for key in bbox:
        if key not in text_kws["bbox"]:
            text_kws["bbox"][key] = bbox[key]
    # plot each text
    text_list = []
    for text, sub_df in data.groupby(anno_col):
        if text_transform is None:
            text = str(text)
        else:
            text = text_transform(text)
        if text.lower() in ["", "nan"]:
            continue
        _x, _y = sub_df[[x, y]].median()

        use_text_kws = copy.deepcopy(text_kws)  # text_kws.copy()
        if isinstance(text_kws["bbox"]["facecolor"], dict):
            use_text_kws["bbox"]["facecolor"] = text_kws["bbox"]["facecolor"].get(
                text, "gray"
            )
        if isinstance(text_kws["color"], dict):
            use_color = text_kws["color"].get(text, "black")
            use_text_kws["color"] = use_color
        if luminance is not None and use_text_kws["bbox"]["facecolor"] is not None:
            lum = _calculate_luminance(use_text_kws["bbox"]["facecolor"])
            if lum > luminance:
                use_text_kws["color"] = "black"
                use_text_kws["bbox"]["edgecolor"] = "black"

        text = ax.text(_x, _y, text, **use_text_kws)
        text_list.append(text)

    if dodge_text:
        try:
            from adjustText import adjust_text

            _dodge_parms = {
                "force_points": (0.02, 0.05),
                "arrowprops": {
                    "arrowstyle": "->",
                    "fc": "black",
                    "ec": "none",
                    "connectionstyle": "angle,angleA=-90,angleB=180,rad=5",
                },
                "autoalign": "xy",
            }
            if dodge_kws is not None:
                _dodge_parms.update(dodge_kws)
            adjust_text(text_list, x=data["x"], y=data["y"], **_dodge_parms)
        except ModuleNotFoundError:
            print(
                "Install adjustText package to dodge text, see its github page for help"
            )
    return


def tight_hue_range(hue_data, portion):
    """Automatic select a SMALLEST data range that covers [portion] of the data."""
    hue_data = hue_data[np.isfinite(hue_data)]
    hue_quantiles = hue_data.quantile(q=np.arange(0, 1, 0.01))
    min_window_right = (
        hue_quantiles.rolling(window=int(portion * 100))
        .apply(lambda i: i.max() - i.min(), raw=True)
        .idxmin()
    )
    min_window_left = max(0, min_window_right - portion)
    vmin, vmax = tuple(hue_data.quantile(q=[min_window_left, min_window_right]))
    if np.isfinite(vmin):
        vmin = max(hue_data.min(), vmin)
    else:
        vmin = hue_data.min()
    if np.isfinite(vmax):
        vmax = min(hue_data.max(), vmax)
    else:
        vmax = hue_data.max()
    return vmin, vmax


def density_contour(
    ax,
    data,
    x,
    y,
    groupby=None,
    c="lightgray",
    single_contour_pad=1,
    linewidth=1,
    linestyles="densely dashed",
    palette=None,
    smooth_sigma=None,
    grid_size=200,
    threshold=None,
):
    """
    Draw density-based contour outlines around groups of points.

    Uses a Gaussian-smoothed 2D histogram to produce clean, smooth outlines.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    data : pd.DataFrame or anndata.AnnData
        Data containing coordinates.
    x, y : str
        Column names for x and y coordinates.
    groupby : str or array-like, optional
        Column name or series defining groups. Each group gets its own contour.
    c : str, default "lightgray"
        Default contour color when palette is not provided.
    single_contour_pad : float, default 1
        Multiplier for auto-computed smooth_sigma. Larger values produce
        contours with more padding around points. Only used when
        smooth_sigma is None.
    linewidth : float, default 1
        Width of the contour lines.
    linestyles : str or tuple, default "densely dashed"
        Linestyle for contour lines. Accepts standard matplotlib names
        ("solid", "dashed", etc.), compound names ("densely dashed",
        "loosely dotted", etc.), or dash tuples like (0, (5, 1)).
    palette : dict, optional
        Mapping of group names to colors. If None, uses ``c`` for all groups.
    smooth_sigma : float, optional
        Gaussian smoothing sigma applied to the 2D histogram (in grid units).
        Controls how smooth the contour is. Larger values produce rounder,
        more generalized outlines; smaller values follow the point
        distribution more closely. If None, auto-computed from the median
        nearest-neighbor distance of the group's points, scaled to produce
        a smooth contour that closely follows the cluster shape.
    grid_size : int, default 200
        Number of bins along each axis for the 2D histogram. Higher values
        give finer resolution but slower computation.
    threshold : float, optional
        Contour level on the normalized [0, 1] density field. Lower values
        draw contours further from the points (looser fit); higher values
        draw contours closer to the dense core (tighter fit). Typical
        range: 0.01 (very loose) to 0.2 (very tight). If None,
        auto-computed to enclose ~95% of points.
    """
    from scipy.ndimage import gaussian_filter

    _data = data.copy()

    if groupby is not None:
        if isinstance(groupby, str):
            _data["groupby"] = data[groupby]
        else:
            _data["groupby"] = groupby
    else:
        _data["groupby"] = "one group"

    # Convert named compound linestyles to dash tuples for contour compatibility
    _named_linestyles = {  # https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
        "loosely dotted": (0, (1, 10)),
        "dotted": (0, (1, 1)),
        "densely dotted": (0, (1, 1)),
        "long dash with offset": (5, (10, 3)),
        "loosely dashed": (0, (5, 10)),
        "dashed": (0, (5, 5)),
        "densely dashed": (0, (5, 1)),
        "loosely dashdotted": (0, (3, 10, 1, 10)),
        "dashdotted": (0, (3, 5, 1, 5)),
        "densely dashdotted": (0, (3, 1, 1, 1)),
        "dashdotdotted": (0, (3, 5, 1, 5, 1, 5)),
        "loosely dashdotdotted": (0, (3, 10, 1, 10, 1, 10)),
        "densely dashdotdotted": (0, (3, 1, 1, 1, 1, 1)),
    }
    _standard_linestyles = {"solid", "dashed", "dashdot", "dotted"}
    if isinstance(linestyles, str) and linestyles not in _standard_linestyles:
        linestyles = _named_linestyles.get(linestyles, linestyles)
    if isinstance(linestyles, tuple):
        linestyles = [linestyles]

    for group, sub_data in _data[[x, y, "groupby"]].groupby("groupby"):
        coords = sub_data.iloc[:, :2].values
        if len(coords) < 5:
            continue
        if palette is None:
            _color = c
        else:
            _color = palette[group] if group in palette else c

        xmin, xmax = zoom_min_max(coords[:, 0].min(), coords[:, 0].max(), 1.15)
        ymin, ymax = zoom_min_max(coords[:, 1].min(), coords[:, 1].max(), 1.15)
        x_range = xmax - xmin
        y_range = ymax - ymin

        # Auto-compute smooth_sigma from median nearest-neighbor distance
        if smooth_sigma is None:
            from sklearn.neighbors import NearestNeighbors

            k = min(5, len(coords) - 1)
            nn = NearestNeighbors(n_neighbors=k + 1)
            nn.fit(coords)
            dists, _ = nn.kneighbors(coords)
            med_nn = np.median(dists[:, 1:])
            # Convert physical distance to grid units and scale
            grid_cell_x = x_range / grid_size
            grid_cell_y = y_range / grid_size
            grid_cell = (grid_cell_x + grid_cell_y) / 2
            _sigma = (med_nn / grid_cell) * 1.5 * single_contour_pad
            _sigma = max(_sigma, 2.0)  # minimum smoothing
        else:
            _sigma = smooth_sigma

        # Build a 2D histogram as occupancy grid, then smooth it
        hist, xedges, yedges = np.histogram2d(
            coords[:, 0],
            coords[:, 1],
            bins=grid_size,
            range=[[xmin, xmax], [ymin, ymax]],
        )
        # Gaussian smooth to get a clean density field
        z = gaussian_filter(hist.T, sigma=_sigma)
        # Normalize to [0, 1]
        z = z / z.max() if z.max() > 0 else z

        # Auto-compute threshold to enclose ~95% of points
        if threshold is None:
            # Sample density values at each point location
            xi = np.clip(
                ((coords[:, 0] - xmin) / x_range * (grid_size - 1)).astype(int),
                0,
                grid_size - 1,
            )
            yi = np.clip(
                ((coords[:, 1] - ymin) / y_range * (grid_size - 1)).astype(int),
                0,
                grid_size - 1,
            )
            point_densities = z[yi, xi]
            # Use the 5th percentile so 95% of points are inside the contour
            _threshold = max(np.percentile(point_densities, 5) * 0.8, 1e-4)
        else:
            _threshold = threshold

        # Create mesh grid centers
        xc = (xedges[:-1] + xedges[1:]) / 2
        yc = (yedges[:-1] + yedges[1:]) / 2
        xx, yy = np.meshgrid(xc, yc)

        ax.contour(
            xx,
            yy,
            z,
            levels=[_threshold],
            colors=_color,
            linewidths=linewidth,
            linestyles=linestyles,
        )
    return


def plot_color_dict_legend(
    D, ax=None, title=None, color_text=True, kws=None, luminance=0.5
):
    """
    plot legned for color dict

    Parameters
    ----------
    D: a dict, key is categorical variable, values are colors.
    ax: axes to plot the legend.
    title: title of legend.
    color_text: whether to change the color of text based on the color in D.
    label_side: right of left.
    kws: kws passed to plt.legend.

    Returns
    -------
    ax.legend

    """
    import matplotlib.patches as mpatches

    if ax is None:
        ax = plt.gca()
    lgd_kws = kws.copy() if kws is not None else {}  # bbox_to_anchor=(x,-0.05)
    lgd_kws.setdefault("frameon", True)
    lgd_kws.setdefault("ncol", 1)
    lgd_kws["loc"] = "upper left"
    lgd_kws.setdefault("borderpad", 0.1 * mm2inch * 72)  # 0.1mm
    lgd_kws.setdefault("markerscale", 1)
    lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
    lgd_kws.setdefault("handlelength", 1)  # font size, units is points
    lgd_kws.setdefault(
        "borderaxespad", 0.1
    )  # The pad between the axes and legend border, in font-size units.
    lgd_kws.setdefault(
        "handletextpad", 0.4
    )  # The pad between the legend handle and text, in font-size units.
    lgd_kws.setdefault(
        "labelspacing", 0.15
    )  # gap height between two Patches,  0.05*mm2inch*72
    lgd_kws.setdefault("columnspacing", 0.5)
    # lgd_kws["bbox_transform"] = ax.figure.transFigure
    lgd_kws.setdefault("bbox_to_anchor", (1, 1))
    lgd_kws.setdefault("title", title)
    lgd_kws.setdefault("markerfirst", True)
    hands = [
        mpatches.Patch(color=c, label=h) for h, c in D.items()
    ]  # kws:?mpatches.Patch; rasterized=True
    ms = lgd_kws.pop("markersize", 10)  # noqa: F841
    L = ax.legend(handles=hands, **lgd_kws)
    L._legend_box.align = "center"
    L.get_title().set_ha("center")
    if color_text:
        for text in L.get_texts():
            try:
                lum = _calculate_luminance(D[text.get_text()])
                if luminance is None:
                    text_color = "black"
                else:
                    text_color = "black" if lum > luminance else D[text.get_text()]
                text.set_color(text_color)
            except AttributeError:
                pass
    # ax.add_artist(L)
    ax.grid(False)
    return L


def plot_marker_legend(
    color_dict=None,
    ax=None,
    title=None,
    color_text=True,
    marker="o",
    kws=None,
    luminance=0.5,
):
    """
    plot legned for different marker

    Parameters
    ----------
    D: a dict, key is categorical variable, values are marker.
    ax: axes to plot the legend.
    title: title of legend.
    color_text: whether to change the color of text based on the color in D.
    label_side: right of left.
    kws: kws passed to plt.legend.

    Returns
    -------
    ax.legend

    """
    if ax is None:
        ax = plt.gca()

    lgd_kws = kws.copy() if kws is not None else {}  # bbox_to_anchor=(x,-0.05)
    lgd_kws.setdefault("frameon", True)
    lgd_kws.setdefault("ncol", 1)
    lgd_kws["loc"] = "upper left"
    # lgd_kws["bbox_transform"] = ax.figure.transFigure
    lgd_kws.setdefault("borderpad", 0.2 * mm2inch * 72)  # 0.1mm
    # lgd_kws.setdefault("markerscale", 1)
    lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
    lgd_kws.setdefault("handlelength", 1)  # font size, units is points
    lgd_kws.setdefault(
        "borderaxespad", 0.1
    )  # The pad between the axes and legend border, in font-size units.
    lgd_kws.setdefault(
        "handletextpad",
        0.4,  # 0.2 * mm2inch * 72
    )  # The pad between the legend handle and text, in font-size units.
    lgd_kws.setdefault(
        "labelspacing", 0.15
    )  # gap height between two Patches,  0.05*mm2inch*72
    lgd_kws.setdefault("columnspacing", 0.5)
    lgd_kws.setdefault("bbox_to_anchor", (1, 1))
    lgd_kws.setdefault("title", title)
    lgd_kws.setdefault("markerfirst", True)

    ms = lgd_kws.pop("markersize", 10)
    L = [
        Line2D(
            [0],
            [0],
            color=color,
            marker=marker,
            linestyle="None",
            markersize=ms,
            label=_lab,
        )
        for _lab, color in color_dict.items()
    ]
    L = ax.legend(handles=L, **lgd_kws)
    ax.figure.canvas.draw()
    L._legend_box.align = "center"
    L.get_title().set_ha("center")
    if color_text:
        for text in L.get_texts():
            try:
                lum = _calculate_luminance(color_dict[text.get_text()])
                if luminance is None:
                    text_color = "black"
                else:
                    text_color = (
                        "black" if lum > luminance else color_dict[text.get_text()]
                    )
                text.set_color(text_color)
            except AttributeError:
                pass
    # ax.add_artist(lgd)
    ax.grid(False)
    return L


# Custom handler for legend: circle text as marker + label
class TextWithCircleHandler(HandlerBase):
    """
    Custom legend handler that renders a numbered circle as the legend marker.

    Used by :func:`plot_text_legend` to create coded legends where each
    category is represented by a colored circle containing a number.

    Parameters
    ----------
    marker_text : str, default ""
        Text to display inside the circle marker (typically a number code).
    label_text : str, default ""
        Label text for the legend entry (category name).
    text_kws : dict, default {}
        Keyword arguments for the marker text, including ``bbox`` for
        the circle style (boxstyle, facecolor, edgecolor, etc.).
    **kwargs
        Additional keyword arguments passed to ``HandlerBase``.
    """

    def __init__(self, marker_text="", label_text="", text_kws={}, **kwargs):
        """Initialize with marker text, label, and text styling kwargs."""
        HandlerBase.__init__(self, **kwargs)
        self.marker_text = marker_text
        self.text_kws = text_kws

    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        """Create a Text artist with circle bbox as the legend marker."""
        self.text_kws.setdefault("fontsize", fontsize)
        # print(self.text_kws)
        shift = 2 * self.text_kws["fontsize"] * 0.65 / 72 / mm2inch
        circ_text = Text(
            xdescent + legend.borderaxespad + shift,
            height / 2,
            self.marker_text,
            **self.text_kws,
        )
        return [circ_text]


def plot_text_legend(
    color_dict,
    code2label,
    ax=None,
    title=None,
    color_text=True,
    boxstyle="Circle",
    marker_pad=0.1,
    legend_kws=None,
    marker_fontsize=4,
    text_kws=None,
    alpha=0.7,
    luminance=0.5,
):
    """
    Plot a coded legend where each category has a numbered circle marker.

    Each legend entry shows a colored circle with a number code followed
    by the category label text. Useful when categories have integer codes
    displayed on the scatter plot.

    Parameters
    ----------
    color_dict : dict
        Mapping of category names to colors (hex strings or RGB tuples).
    code2label : dict
        Mapping of integer codes to category label strings,
        e.g. ``{1: "Excitatory", 2: "Inhibitory"}``.
    ax : matplotlib.axes.Axes, optional
        Axes to attach the legend to. Uses current axes if None.
    title : str, optional
        Legend title.
    color_text : bool, default True
        If True, color the label text to match the category color
        (for labels with luminance below ``luminance``).
    boxstyle : str, default "Circle"
        Shape of the marker box. Options: "Circle", "Round", "Square".
    marker_pad : float, default 0.1
        Padding inside the marker box (passed to matplotlib FancyBboxPatch).
    legend_kws : dict, optional
        Additional keyword arguments passed to ``ax.legend()``, such as
        ``ncol``, ``fontsize``, ``bbox_to_anchor``, etc.
    marker_fontsize : float, default 4
        Font size of the number text inside the circle marker.
    text_kws : dict, optional
        Keyword arguments for the marker text rendering, including ``bbox``
        sub-dict for circle styling (edgecolor, linewidth, facecolor, etc.).
    alpha : float, default 0.7
        Opacity of the circle marker background.
    luminance : float, default 0.5
        Luminance threshold for deciding text color. Labels with color
        luminance above this value are shown in black; darker colors
        are shown in their own color.
    """
    lgd_kws = (
        legend_kws.copy() if legend_kws is not None else {}
    )  # bbox_to_anchor=(x,-0.05)
    lgd_kws.setdefault("frameon", True)
    lgd_kws.setdefault("ncol", 1)
    lgd_kws["loc"] = "upper left"
    # lgd_kws["bbox_transform"] = ax.figure.transFigure
    lgd_kws.setdefault("borderpad", 0.1 * mm2inch * 72)  # 0.1mm
    # lgd_kws.setdefault("markerscale", 1)
    lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
    lgd_kws.setdefault("handlelength", 1)  # font size, units is points
    lgd_kws.setdefault(
        "borderaxespad", 0.5
    )  # The pad between the axes and legend border, in font-size units.
    lgd_kws.setdefault(
        "handletextpad", 0.4
    )  # The pad between the legend handle and text, in font-size units.
    lgd_kws.setdefault(
        "labelspacing", 0.2
    )  # gap height between two row of legend,  0.05*mm2inch*72
    lgd_kws.setdefault("columnspacing", 0.5)
    lgd_kws.setdefault("bbox_to_anchor", (1, 1))
    lgd_kws.setdefault("title", title)
    lgd_kws.setdefault("markerfirst", True)
    ms = lgd_kws.pop("markersize", 10)  # noqa:F841

    # text_kws
    if text_kws is None:
        text_kws = {}
    default_marker_text_kws = dict(
        bbox=dict(
            boxstyle=f"{boxstyle},pad={marker_pad}",  # Square, Circle, Round
            edgecolor="black",
            linewidth=0.4,
            fill=True,
            facecolor="white",
            alpha=alpha,
        ),
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=marker_fontsize,
        color="black",
    )
    for k in default_marker_text_kws:
        if k == "bbox":
            if k not in text_kws:
                text_kws["bbox"] = default_marker_text_kws[k]
            else:
                for k1 in default_marker_text_kws["bbox"].keys():
                    if k1 not in text_kws["bbox"]:
                        text_kws["bbox"].setdefault(
                            k1, default_marker_text_kws["bbox"][k1]
                        )
        if k not in text_kws:
            text_kws.setdefault(k, default_marker_text_kws[k])

    # Create handles and handlers
    handles = []
    handler_map = {}
    for code in sorted([int(i) for i in code2label.keys()]):
        label = code2label[code]
        code_text = str(code)
        handle = Line2D([], [], linestyle=None, label=label)
        handles.append(handle)
        color = color_dict.get(label, "black")
        text_kws1 = copy.deepcopy(text_kws)
        text_kws1["bbox"]["facecolor"] = color
        lum = _calculate_luminance(color)
        if lum <= 0.1:  # for black-like color, use white marker text
            text_kws1["color"] = "white"
        # print(bbox)
        handler_map[handle] = TextWithCircleHandler(
            marker_text=code_text,
            label_text=label,
            text_kws=text_kws1,
        )

    # Draw custom legend
    L = ax.legend(handles=handles, handler_map=handler_map, **lgd_kws)
    L._legend_box.align = "center"
    L.get_title().set_ha("center")
    ax.figure.canvas.draw()
    if color_text:
        for text in L.get_texts():
            try:
                lum = _calculate_luminance(color_dict[text.get_text()])
                if luminance is None:
                    text_color = "black"
                else:
                    text_color = (
                        "black" if lum > luminance else color_dict[text.get_text()]
                    )
                text.set_color(text_color)
            except AttributeError:
                pass
    # ax.add_artist(lgd)
    # ax.grid(False)


def plot_cmap_legend(
    cax=None,
    ax=None,
    cmap="turbo",
    label=None,
    kws=None,
    labelsize=6,
    linewidth=0.5,
    ticklabel_size=4,
    ticklabel_pad=1,
):
    """
    Plot legend for cmap.

    Parameters
    ----------
    cax : Axes into which the colorbar will be drawn.
    ax :  axes to anchor.
    cmap : turbo, hsv, Set1, Dark2, Paired, Accent,tab20,exp1,exp2,meth1,meth2
    label : title for legend.
    kws : dict
            kws passed to plt.colorbar (matplotlib.figure.Figure.colorbar).

    Returns
    -------
    cbar: axes of legend

    """
    label = "" if label is None else label
    cbar_kws = {} if kws is None else kws.copy()
    cbar_kws.setdefault("label", label)
    # cbar_kws.setdefault("aspect",3)
    cbar_kws.setdefault("orientation", "vertical")
    # cbar_kws.setdefault("use_gridspec", True)
    # cbar_kws.setdefault("location", "bottom")
    cbar_kws.setdefault("fraction", 1)
    cbar_kws.setdefault("shrink", 1)
    cbar_kws.setdefault("pad", 0)
    cbar_kws.setdefault("extend", "both")
    cbar_kws.setdefault("extendfrac", 0.1)
    # print(cbar_kws,kws)
    # print(type(cax))
    vmax = cbar_kws.pop("vmax", 1)
    vmin = cbar_kws.pop("vmin", 0)
    vcenter = (vmax + vmin) / 2
    center = cbar_kws.pop("center", None)
    if center is None:
        center = vcenter
        m = plt.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap
        )
    else:
        m = plt.cm.ScalarMappable(
            norm=matplotlib.colors.TwoSlopeNorm(center, vmin=vmin, vmax=vmax), cmap=cmap
        )
    cbar_kws.setdefault("ticks", [vmin, center, vmax])
    cax.yaxis.set_label_position("right")
    cax.yaxis.set_ticks_position("right")
    cbar = ax.figure.colorbar(m, cax=cax, **cbar_kws)  # use_gridspec=True
    # Fix while-line bug
    cbar.solids.set_visible(False)
    gradient = np.linspace(vmin, vmax, 1024).reshape(-1, 1)
    cax.imshow(
        gradient,
        aspect="auto",
        cmap=m.cmap,
        norm=m.norm,
        origin="lower",
        extent=(0, 1, vmin, vmax),
        interpolation="bilinear",
        zorder=cbar.solids.get_zorder(),
    )
    cax.set_xlim(0, 1)
    cax.set_ylim(vmin, vmax)
    cbar.ax.tick_params(
        labelsize=ticklabel_size,
        size=ticklabel_size,
        width=linewidth,
        pad=ticklabel_pad,
    )  # size is for ticks, labelsize is for the number on ticks (ticklabels), pad is gap between tick and label
    cbar.ax.yaxis.label.set_fontsize(labelsize)  # colorbar title fontsize
    cbar.ax.grid(False)
    return cbar


def normalize_mc_by_cell(
    use_adata, normalize_per_cell=True, clip_norm_value=10, verbose=1, hypo_score=False
):
    """
    Normalize methylation fraction per cell by prior mean.

    Divides each cell's methylation values by its prior mean (stored in
    ``obs['prior_mean']``), which is typically determined by alpha and
    beta parameters of a Beta distribution.

    Parameters
    ----------
    use_adata : anndata.AnnData
        AnnData object with methylation fractions in ``.X`` and
        ``prior_mean`` column in ``.obs``.
    normalize_per_cell : bool, default True
        Whether to perform normalization. If False, returns the input
        unchanged.
    clip_norm_value : float or None, default 10
        Maximum value to clip normalized values to. Set to None to
        disable clipping.
    verbose : int, default 1
        Verbosity level. 0 = silent, 1 = print status messages.
    hypo_score : bool, default False
        If True, compute hypomethylation score (prior_mean / value)
        instead of the standard normalized fraction (value / prior_mean).
        Higher hypo_score means more hypomethylated.

    Returns
    -------
    anndata.AnnData
        The input AnnData with ``.X`` modified in-place.
        Sets ``uns['normalize_per_cell'] = True`` after normalization.
    """
    from scipy.sparse import issparse

    normalized_flag = use_adata.uns.get("normalize_per_cell", False)
    if (
        normalize_per_cell and not normalized_flag
    ):  # divide frac by prior mean (determined by alpha and beta) for each cell
        # get normalized X
        cols = use_adata.obs.columns.tolist()
        if "prior_mean" in cols:
            if verbose > 0:
                print("Normalizing cell level fraction by alpha and beta (prior_mean)")
            na_sum = use_adata.to_df().isna().sum().sum()
            if na_sum > 0:
                D = use_adata.obs.prior_mean.to_dict()
                if not hypo_score:
                    use_adata.X = (
                        use_adata.to_df()
                        .apply(lambda x: x.fillna(D[x.name]) / D[x.name], axis=1)
                        .values
                    )
                else:  # hypo score, the larger the value, the more hypomethylated
                    use_adata.X = (
                        use_adata.to_df()
                        .apply(lambda x: D[x.name] / x.fillna(D[x.name]), axis=1)
                        .values
                    )
            else:
                if (
                    not hypo_score
                ):  # the smaller the value, the lower of the methylation fraction
                    use_adata.X = (
                        use_adata.X / use_adata.obs.prior_mean.values[:, None]
                    )  # range = [0,1,10]
                else:  # hypo-score
                    use_adata.X = use_adata.obs.prior_mean.values[:, None] / use_adata.X
            if clip_norm_value is not None:
                if issparse(use_adata.X):
                    X = use_adata.X.toarray()
                else:
                    X = use_adata.X
                use_adata.X = np.clip(X, None, clip_norm_value)
            use_adata.uns["normalize_per_cell"] = True
        else:
            if verbose > 0:
                print("'prior_mean' not found in obs")
    elif normalize_per_cell and normalized_flag:
        # logger.info("Input adata is already normalized, skip normalize_per_cell !")
        print(
            "Input adata is already normalized, skip normalize_per_cell !"
        )  # TODO: use logger
    else:
        pass
    return use_adata


def parse_gtf(gtf="gencode.v43.annotation.gtf", outfile=None):
    """
    Parse a GENCODE/Ensembl GTF annotation file into a DataFrame.

    Extracts gene-level information including gene_id, gene_name,
    gene_type, genomic coordinates, and strand.

    Parameters
    ----------
    gtf : str, default "gencode.v43.annotation.gtf"
        Path to the GTF file. Supports ``~`` expansion.
    outfile : str or None, default None
        If provided, saves the parsed result as a tab-separated file
        at this path. If None, returns the DataFrame directly.

    Returns
    -------
    pd.DataFrame or None
        DataFrame with columns: chrom, beg, end, gene_name, gene_id,
        strand, gene_type. Returned only when ``outfile`` is None.
    """
    df = pd.read_csv(
        os.path.expanduser(gtf),
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 2, 3, 4, 6, 8],
        names=["chrom", "record_type", "beg", "end", "strand", "information"],
    )
    cols = ["gene_id", "gene_type", "gene_name"]

    def parse_info(x):
        """Parse a GTF info field string into a key-value dict."""
        x = x.replace('"', "")
        D = {}
        for item in x.strip().rstrip(";").split(";"):
            k, v = item.strip().split(" ")
            D[k.strip()] = v.strip()
        return D

    df["info_dict"] = df.information.apply(parse_info)
    for col in cols:
        df[col] = df.info_dict.apply(lambda x: x.get(col, ""))
    df = df.loc[
        :, ["chrom", "beg", "end", "gene_name", "gene_id", "strand", "gene_type"]
    ].drop_duplicates()
    if outfile is None:
        return df  # 'chrom','beg','end','gene_name','gene_id','strand','gene_type'
    else:
        df.to_csv(os.path.expanduser(outfile), sep="\t", index=False)


def categorical_scatter(
    data,
    ax=None,
    basis="umap",
    x=None,
    y=None,  # coords
    hue=None,
    palette="auto",
    color=None,  # color
    text_anno=None,
    text_kws=None,
    luminance=None,
    text_transform=None,
    dodge_text=False,
    dodge_kws=None,  # text annotation
    show_legend=False,
    legend_kws=None,  # legend
    s="auto",
    size=None,
    sizes=None,  # sizes is a dict
    size_norm=None,
    size_portion=0.95,
    axis_format="tiny",
    max_points=50000,
    labelsize=4,
    linewidth=0.5,
    zoomxy=1.05,
    outline=None,
    outline_pad=3,
    alpha=0.7,
    outline_kws=None,
    scatter_kws=None,
    rasterized="auto",
    coding=False,
    coded_marker=True,
    legend_color_text=True,
    rectangle_marker=False,
    marker_fontsize=4,
    marker_pad=0.1,
):
    """
    This function was copied from ALLCools and made some modifications.
    Plot categorical scatter plot with versatile options.

    Parameters
    ----------
    rasterized
            Whether to rasterize the figure.
    return_fig
            Whether to return the figure.
    size_portion
            The portion of the figure to be used for the size norm.
    data
            Dataframe that contains coordinates and categorical variables
    ax
            this function do not generate ax, must provide an ax
    basis
            coords name, if provided, will automatically search for x and y
    x
            x coord name
    y
            y coord name
    hue : str
            categorical col name or series for color hue.
    palette : str or dict
            palette for color hue.
    color
            specify single color for all the dots
    text_anno
            categorical col name or series for text annotation.
    text_kws
            kwargs pass to plt.text, see: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
            including bbox, to see parameter for bbox, go to: https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.FancyBboxPatch.html#matplotlib.patches.FancyBboxPatch
            commonly used parameters are::

                    text_kws=dict(fontsize=5,fontweight='black',
                                            color='black', # color could be a dict, keys are text to be annotated
                                            bbox=dict(boxstyle='round',edgecolor=(0.5, 0.5, 0.5, 0.2),fill=False,
                                                                    facecolor=(0.8, 0.8, 0.8, 0.2), # facecolor could also be a dict
                                                                    alpha=1,linewidth=0.5)
                                            )
    text_transform
            transform for text annotation.
    dodge_text
            whether to dodge text annotation.
    dodge_kws
            kwargs for dodge text annotation.
    show_legend
            whether to show legend.
    legend_kws
            kwargs for legend.
    s
            single size value of all the dots.
    size
            mappable size of the dots.
    sizes
            mapping size to the sizes value.
    size_norm
            normalize size range for mapping.
    axis_format
            axis format.
    max_points
            maximum number of points to plot.
    labelsize
            label size pass to `ax.text`
    linewidth
            line width pass to `ax.scatter`
    zoomxy
            zoom factor for x and y-axis.
    outline
            categorical col name or series for outline.
    outline_pad
            outline padding.
    outline_kws
            kwargs for outline.
    scatter_kws
            kwargs for scatter.

    Returns
    -------
    if return_fig is True, return the figure and axes.
    else, return None.
    """
    if ax is None:
        # fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
        ax = plt.gca()
    # add coords
    _data, x, y = _extract_coords(data, basis, x, y)
    # _data has 2 cols: "x" and "y", index are obs_names

    # down sample plot data if needed.
    if max_points is not None:
        if _data.shape[0] > max_points:
            _data = _density_based_sample(
                _data, seed=1, size=max_points, coords=["x", "y"]
            )
    n_dots = _data.shape[0]

    # determine rasterized
    if rasterized == "auto":
        if n_dots > 200:
            rasterized = True
        else:
            rasterized = False

    # auto size if user didn't provide one
    if s == "auto":
        s = _auto_size(ax, n_dots)

    # default scatter options
    _scatter_kws = {
        "linewidth": 0,
        "s": s,
        "legend": None,
        "palette": palette,
        "rasterized": rasterized,
    }
    if color is not None:
        if hue is not None:
            raise ValueError("Only one of color and hue can be provided")
        _scatter_kws["color"] = color
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)

    # deal with color
    palette_dict = None
    if hue is not None:
        if isinstance(hue, str):
            _data["hue"] = _take_data_series(data, hue)
        else:
            _data["hue"] = hue.copy()
        _data["hue"] = _data["hue"].astype("category").cat.remove_unused_categories()

        # if the object has get_palette method, use it (AnnotZarr)
        palette = _scatter_kws["palette"]
        # deal with other color palette
        if palette_dict is None:
            if isinstance(palette, str) or isinstance(palette, list):
                palette_dict = level_one_palette(
                    _data["hue"], order=None, palette=palette
                )
            elif isinstance(palette, dict):
                palette_dict = palette
            else:
                raise TypeError(
                    f"Palette can only be str, list or dict, got {type(palette)}"
                )
        _scatter_kws["palette"] = palette_dict
    else:
        _scatter_kws["palette"] = None

    # deal with size
    if size is not None:
        if isinstance(size, str):
            _data["size"] = _take_data_series(data, size).astype(float)
        else:
            _data["size"] = size.astype(float)
        size = "size"

        if size_norm is None:
            # get the smallest range that include "size_portion" of data
            size_norm = tight_hue_range(_data["size"], size_portion)

            # snorm is the normalizer for size
            size_norm = matplotlib.colors.Normalize(
                vmin=size_norm[0], vmax=size_norm[1]
            )

        # discard s from _scatter_kws and use size in sns.scatterplot
        s = _scatter_kws.pop("s")
        if sizes is None:
            sizes = (min(s, 1), s)

    sns.scatterplot(
        x="x",
        y="y",
        data=_data,
        ax=ax,
        hue="hue" if hue is not None else None,
        size=size,
        sizes=sizes,
        size_norm=size_norm,
        **_scatter_kws,
    )

    # deal with text annotation
    code2label = None
    if text_anno is not None:
        # data
        if isinstance(text_anno, str):
            _data["text_anno"] = _take_data_series(data, text_anno)
        else:
            _data["text_anno"] = text_anno.copy()
        if str(_data["text_anno"].dtype) == "category":
            _data["text_anno"] = _data["text_anno"].cat.remove_unused_categories()

        # text kws
        text_kws = {} if text_kws is None else text_kws
        default_text_kws = dict(
            color="white",  # color for the text, could be a dict, keys are text to be annotated
            fontweight="bold",  # fontsize=labelsize,
            bbox=dict(
                facecolor=palette_dict,  # if None, use default color
                boxstyle="round",  # ellipse, round
                edgecolor="white",
                fill=True,
                linewidth=linewidth,
                alpha=alpha,
            ),
        )
        # coding & coded_marker
        text_anno = "text_anno"
        if (coding is not None) and (coding != False):  # noqa: E712
            if coding == True:  # noqa: E712
                _data["code"] = _data[
                    "hue"
                ].cat.codes  # int; hue already cleaned, codes are contiguous
            else:
                assert isinstance(coding, str)
                _data["code"] = _take_data_series(data, coding)
                _data = _data.loc[_data["code"].notna()]
                _data["code"] = _data["code"].astype(int)
            text_anno = "code"
            _data["color"] = [
                palette_dict.get(v) for v in _data["hue"]
            ]  # avoid Categorical.map() pandas 3 bug
            code2label = (
                _data.loc[:, ["code", "hue"]]
                .drop_duplicates()
                .set_index("code")
                .hue.to_dict()
            )
            _data["code"] = _data["code"].astype(str)
            code_colors = (
                _data.loc[:, ["code", "color"]]
                .drop_duplicates()
                .set_index("code")
                .color.to_dict()
            )
            default_text_kws["bbox"]["facecolor"] = (
                code_colors  # background colors for text annotation
            )
            default_text_kws["bbox"]["boxstyle"] = "circle"
        for k in default_text_kws:
            if k != "bbox":
                text_kws.setdefault(k, default_text_kws[k])
            else:
                if "bbox" not in text_kws:
                    text_kws["bbox"] = {}
                for k1 in default_text_kws["bbox"]:
                    text_kws["bbox"].setdefault(k1, default_text_kws["bbox"][k1])

        _text_anno_scatter(
            data=_data[["x", "y", text_anno]],
            ax=ax,
            x="x",
            y="y",
            dodge_text=dodge_text,
            dodge_kws=dodge_kws,
            text_transform=text_transform,
            anno_col=text_anno,
            text_kws=text_kws,
            luminance=luminance,
        )

    # deal with outline
    if outline is not None:
        if isinstance(outline, str):
            _data["outline"] = _take_data_series(data, outline)
        else:
            _data["outline"] = outline.copy()
        _outline_kws = {
            "linewidth": linewidth,
            "palette": None,
            "c": "lightgray",
            "single_contour_pad": outline_pad,
        }
        if outline_kws is not None:
            _outline_kws.update(outline_kws)
        density_contour(
            ax=ax, data=_data, x="x", y="y", groupby="outline", **_outline_kws
        )

    # clean axis
    if axis_format == "tiny":
        _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
    elif (axis_format == "empty") or (axis_format is None):
        despine(ax=ax, left=True, bottom=True)
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    else:
        pass

    # deal with legend
    if show_legend and (hue is not None):
        n_hue = len(palette_dict)
        ncol = 1 if n_hue <= 40 else 2 if n_hue <= 100 else 3
        if legend_kws is None:
            legend_kws = {}
        default_lgd_kws = dict(
            ncol=ncol,
            fontsize=labelsize,
            bbox_to_anchor=(1, 1),
            loc="upper left",
            # borderpad=0.4, # pad between marker (text) and border
            # labelspacing=0.2, #The vertical space between the legend entries, in font-size units
            # handleheight=0.5, #The height of the legend handles, in font-size units.
            # handletextpad=0.2, # The pad between the legend handle (marker) and text, in font-size units.
            # borderaxespad=0.3, # The pad between the Axes and legend border, in font-size units
            # columnspacing=0.2, #The spacing between columns, in font-size units
            markersize=labelsize,  # legend_kws["fontsize"],
        )
        for k in default_lgd_kws:
            legend_kws.setdefault(k, default_lgd_kws[k])

        exist_hues = _data["hue"].unique()
        color_dict = {
            hue_name: color
            for hue_name, color in palette_dict.items()
            if hue_name in exist_hues
        }

        if code2label is not None and coded_marker:
            boxstyle = "Circle" if not rectangle_marker else "Round"
            plot_text_legend(
                color_dict,
                code2label,
                ax,
                title=hue,
                color_text=legend_color_text,
                boxstyle=boxstyle,
                marker_pad=marker_pad,
                legend_kws=legend_kws,
                marker_fontsize=marker_fontsize,
                alpha=alpha,
                luminance=luminance,
            )
        else:
            if rectangle_marker:
                ## plot Patch legend (rectangle marker)
                plot_color_dict_legend(
                    D=color_dict,
                    ax=ax,
                    title=hue,
                    color_text=legend_color_text,
                    kws=legend_kws,
                    luminance=luminance,
                )
            else:
                # plot marker legend (for example, circle marker)
                plot_marker_legend(
                    color_dict=color_dict,
                    ax=ax,
                    title=hue,
                    color_text=legend_color_text,
                    marker="o",
                    kws=legend_kws,
                    luminance=luminance,
                )

    if zoomxy is not None:
        zoom_ax(ax, zoomxy)

    return _data


def get_cmap(cmap):
    """
    Get a matplotlib colormap by name, including custom ones.

    .. deprecated::
        Use :func:`adataviz.palettes.get_cmap` instead.

    Parameters
    ----------
    cmap : str
        Name of the colormap. Supports standard matplotlib names
        ("viridis", "turbo", etc.) and custom adataviz colormaps
        ("parula", "exp1", "exp2", "meth1", "meth2").

    Returns
    -------
    matplotlib.colors.Colormap
        The requested colormap object.
    """
    from . import palettes as _palettes

    warnings.warn(
        "adataviz.utils.get_cmap is deprecated; use adataviz.palettes.get_cmap instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _palettes.get_cmap(cmap)


def continuous_scatter(
    data,
    ax=None,
    basis="umap",
    x=None,
    y=None,
    scatter_kws=None,
    hue=None,
    hue_norm=None,
    hue_portion=0.95,
    color=None,
    cmap="viridis",
    colorbar=True,
    size=None,
    size_norm=None,
    size_portion=0.95,
    sizes=None,
    sizebar=True,
    text_anno=None,
    dodge_text=False,
    dodge_kws=None,
    text_kws=None,
    luminance=0.48,
    text_transform=None,
    axis_format="tiny",
    max_points=50000,
    s="auto",
    labelsize=6,
    ticklabel_size=4,
    linewidth=0.5,
    zoomxy=1.05,
    outline=None,
    outline_kws=None,
    outline_pad=2,
    return_fig=False,
    rasterized="auto",
    cbar_kws=None,
    cbar_width=3,
):
    """
    Plot scatter on given adata.

    Parameters
    ----------
    data : _type_
            _description_
    ax : _type_, optional
            _description_, by default None
    basis : str, optional
            _description_, by default "umap"
    x : _type_, optional
            _description_, by default None
    y : _type_, optional
            _description_, by default None
    scatter_kws : _type_, optional
            _description_, by default None
    hue : _type_, optional
            _description_, by default None
    hue_norm : _type_, optional
            _description_, by default None
    hue_portion : float, optional
            _description_, by default 0.95
    color : _type_, optional
            _description_, by default None
    cmap : str, optional
            _description_, by default "viridis"
    colorbar : bool, optional
            _description_, by default True
    size : _type_, optional
            _description_, by default None
    size_norm : _type_, optional
            _description_, by default None
    size_portion : float, optional
            _description_, by default 0.95
    sizes : _type_, optional
            _description_, by default None
    sizebar : bool, optional
            _description_, by default True
    text_anno : _type_, optional
            _description_, by default None
    dodge_text : bool, optional
            _description_, by default False
    dodge_kws : _type_, optional
            _description_, by default None
    text_kws : _type_, optional
            _description_, by default None
    luminance : float, optional
            _description_, by default 0.48
    text_transform : _type_, optional
            _description_, by default None
    axis_format : str, optional
            _description_, by default "tiny"
    max_points : int, optional
            _description_, by default 50000
    s : str, optional
            _description_, by default "auto"
    labelsize : int, optional
            _description_, by default 6
    ticklabel_size : int, optional
            _description_, by default 4
    linewidth : float, optional
            _description_, by default 0.5
    zoomxy : float, optional
            _description_, by default 1.05
    outline : _type_, optional
            _description_, by default None
    outline_kws : _type_, optional
            _description_, by default None
    outline_pad : int, optional
            _description_, by default 2
    return_fig : bool, optional
            _description_, by default False
    rasterized : str, optional
            _description_, by default "auto"
    cbar_kws : _type_, optional
            _description_, by default None
    cbar_width : int, optional
            width of colorbar, by default 3 mm

    Returns
    -------
    _type_
            _description_

    Raises
    ------
    ValueError
            _description_
    TypeError
            _description_
    """
    import seaborn as sns
    from matplotlib.cm import ScalarMappable

    # init figure if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    else:
        fig = None

    # add coords
    _data, x, y = _extract_coords(data, basis, x, y)
    # _data has 2 cols: "x" and "y"

    # down sample plot data if needed.
    if max_points is not None:
        if _data.shape[0] > max_points:
            _data = _density_based_sample(
                _data, seed=1, size=max_points, coords=["x", "y"]
            )
    n_dots = _data.shape[0]

    # determine rasterized
    if rasterized == "auto":
        if n_dots > 200:
            rasterized = True
        else:
            rasterized = False

    # auto size if user didn't provide one
    if s == "auto":
        s = _auto_size(ax, n_dots)

    # default scatter options
    _scatter_kws = {"linewidth": 0, "s": s, "legend": None, "rasterized": rasterized}
    if color is not None:
        if hue is not None:
            raise ValueError("Only one of color and hue can be provided")
        _scatter_kws["color"] = color
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)

    # deal with color
    if hue is not None:
        if isinstance(hue, str):
            _data["hue"] = _take_data_series(data, hue).astype(float)
            # colorbar_label = hue
        else:
            _data["hue"] = hue.astype(float)
            # colorbar_label = hue.name

        if hue_norm is None:
            # get the smallest range that include "hue_portion" of data
            # hue_norm = tight_hue_range(_data["hue"], hue_portion)
            hue_norm = (
                _data["hue"].quantile(1 - hue_portion),
                _data["hue"].quantile(hue_portion),
            )
        # cnorm is the normalizer for color
        cnorm = matplotlib.colors.Normalize(vmin=hue_norm[0], vmax=hue_norm[1])
        if isinstance(cmap, str):
            # from here, cmap become colormap object
            cmap = copy.copy(get_cmap(cmap))
            cmap.set_bad(color=(0.5, 0.5, 0.5, 0.5))
        else:
            if not isinstance(cmap, ScalarMappable):
                raise TypeError(
                    f"cmap can only be str or ScalarMappable, got {type(cmap)}"
                )
    else:
        hue_norm = None
        cnorm = None
        # colorbar_label = "" #TODO: use colorbar label downstream

    # deal with size
    if size is not None:
        if isinstance(size, str):
            _data["size"] = _take_data_series(data, size).astype(float)
        else:
            _data["size"] = size.astype(float)
        size = "size"

        if size_norm is None:
            # get the smallest range that include "size_portion" of data
            size_norm = tight_hue_range(_data["size"], size_portion)

            # snorm is the normalizer for size
            size_norm = matplotlib.colors.Normalize(
                vmin=size_norm[0], vmax=size_norm[1]
            )

        # replace s with sizes
        s = _scatter_kws.pop("s")
        if sizes is None:
            sizes = (min(s, 1), s)
    else:
        size_norm = None
        sizes = None

    sns.scatterplot(
        x="x",
        y="y",
        data=_data,
        hue="hue" if hue is not None else None,
        palette=cmap if hue is not None else None,
        hue_norm=cnorm,
        size=size,
        sizes=sizes,
        size_norm=size_norm,
        ax=ax,
        **_scatter_kws,
    )

    if text_anno is not None:
        if isinstance(text_anno, str):
            _data["text_anno"] = _take_data_series(data, text_anno)
        else:
            _data["text_anno"] = text_anno
        if str(_data["text_anno"].dtype) == "category":
            _data["text_anno"] = _data["text_anno"].cat.remove_unused_categories()

        _text_anno_scatter(
            data=_data[["x", "y", "text_anno"]],
            ax=ax,
            x="x",
            y="y",
            dodge_text=dodge_text,
            dodge_kws=dodge_kws,
            text_transform=text_transform,
            anno_col="text_anno",
            text_kws=text_kws,
            luminance=luminance,
        )

    # deal with outline
    if outline:
        if isinstance(outline, str):
            _data["outline"] = _take_data_series(data, outline)
        else:
            _data["outline"] = outline
        _outline_kws = {
            "linewidth": linewidth,
            "palette": None,
            "c": "lightgray",
            "single_contour_pad": outline_pad,
        }
        if outline_kws is not None:
            _outline_kws.update(outline_kws)
        density_contour(
            ax=ax, data=_data, x="x", y="y", groupby="outline", **_outline_kws
        )

    # clean axis
    if axis_format == "tiny":
        _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
    elif (axis_format == "empty") or (axis_format is None):
        despine(ax=ax, left=True, bottom=True)
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    else:
        pass

    return_axes = [ax]

    # make color bar
    if colorbar and (hue is not None):
        # small ax for colorbar
        # default_cbar_kws=dict(loc="upper left", borderpad=0,width="3%", height="20%") #bbox_to_anchor=(1,1)
        if cbar_kws is None:
            cbar_kws = {}
        # for k in default_cbar_kws:
        #     if k not in cbar_kws:
        #         cbar_kws[k]=default_cbar_kws[k]

        mm2inch = 1 / 25.4
        space = 0
        legend_width = (
            cbar_width * mm2inch * ax.figure.dpi / ax.figure.get_window_extent().width
        )  # mm to px to fraction
        pad = (
            space + ax.yaxis.labelpad * 1.2 * ax.figure.dpi / 72
        ) / ax.figure.get_window_extent().width
        # labelpad unit is points
        left = ax.get_position().x1 + pad
        ax_legend = ax.figure.add_axes(
            [
                left,
                ax.get_position().height * 0.8,
                legend_width,
                ax.get_position().height * 0.2,
            ]
        )  # left, bottom, width, height
        # print("test:",hue_norm)
        # cbar_kws.setdefault('vmin',hue_norm[0])
        # cbar_kws.setdefault('vmax',hue_norm[1])
        cbar_kws["vmin"] = hue_norm[0]
        cbar_kws["vmax"] = hue_norm[1]
        _ticklabel_pad = cbar_kws.pop("ticklabel_pad", 1)
        cbar = plot_cmap_legend(
            ax=ax,
            cax=ax_legend,
            cmap=cmap,
            label=hue,
            kws=cbar_kws.copy(),
            labelsize=labelsize,
            linewidth=linewidth,
            ticklabel_size=ticklabel_size,
            ticklabel_pad=_ticklabel_pad,
        )
        return_axes.append([ax_legend, cbar])

    # make size bar
    if sizebar and (size is not None):
        # TODO plot dot size bar
        pass

    if zoomxy is not None:
        zoom_ax(ax, zoomxy)

    if return_fig:
        return (fig, tuple(return_axes)), _data
    else:
        return
