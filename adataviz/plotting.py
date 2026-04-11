import os
import pandas as pd
import anndata
import matplotlib.pylab as plt
import numpy as np
from .utils import normalize_mc_by_cell, categorical_scatter, continuous_scatter
from .tools import load_adata, load_obs, load_color_palette
from loguru import logger as logger
# logger.add(sys.stderr, level="DEBUG")
# logger.add(sys.stderr, level="ERROR")


def use_scientific_style():
    """Configure matplotlib rcParams for publication-quality figures.

    Sets font to Arial (8–11 pt), high DPI (300), tight/transparent saves,
    and Type 42 fonts for PDF/PS compatibility. Call once before plotting.
    """
    import matplotlib.pylab as plt

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": "Arial",
            # Base font size (set so text is ~7–8 pt at final print size)
            "font.size": 8,  # main text (labels, ticks, legend)
            "axes.labelsize": 9,  # axis labels (x/y)
            "axes.titlesize": 10,  # panel titles / figure titles
            "figure.titlesize": 11,
            "xtick.labelsize": 8,  # x-tick labels
            "ytick.labelsize": 8,  # y-tick labels
            "legend.fontsize": 8,  # legend text
            "legend.title_fontsize": 9,  # legend title (if used)
            # Lines and elements for clarity
            "lines.linewidth": 1.2,  # data lines
            "axes.linewidth": 1.0,  # axis spines
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            # General figure appearance
            "figure.dpi": 300,  # high resolution for export
            "savefig.dpi": 300,  # when using plt.savefig()
            "figure.figsize": (
                6.5,
                4.5,
            ),  # example starting size (adjust to your needs; e.g., ~17 cm wide for full page)
            # 'figure.constrained_layout.use': True,
            "figure.autolayout": True,
            "savefig.transparent": True,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    # plt.rcParams.keys()


def _assign_default_colors(adata, key):
    """Assign default colors for a categorical obs column, mimicking scanpy's logic."""
    from matplotlib.colors import to_hex
    from matplotlib import rcParams

    categories = adata.obs[key].cat.categories
    n = len(categories)

    # scanpy's default palettes (from scanpy.plotting.palettes)
    _default_20 = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
        "#aec7e8",
        "#ffbb78",
        "#98df8a",
        "#ff9896",
        "#c5b0d5",
        "#c49c94",
        "#f7b6d2",
        "#c7c7c7",
        "#dbdb8d",
        "#9edae5",
    ]
    _default_28 = _default_20 + [
        "#393b79",
        "#637939",
        "#8c6d31",
        "#843c39",
        "#7b4173",
        "#5254a3",
        "#bd9e39",
        "#ad494a",
    ]
    _default_102 = [
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
        "#6A3A4C",
        "#83AB58",
        "#001C1E",
        "#D1F7CE",
        "#004B28",
        "#C8D0F6",
    ]

    # pick palette based on number of categories (same as scanpy)
    default_cycle = rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(default_cycle) >= n:
        palette = [to_hex(c) for c in default_cycle[:n]]
    elif n <= 20:
        palette = [to_hex(c) for c in _default_20[:n]]
    elif n <= 28:
        palette = [to_hex(c) for c in _default_28[:n]]
    elif n <= len(_default_102):
        palette = [to_hex(c) for c in _default_102[:n]]
    else:
        palette = ["grey"] * n

    adata.uns[f"{key}_colors"] = palette
    return {cat: col for cat, col in zip(categories, palette)}


def interactive_embedding(
    adata=None,
    obs=None,
    variable=None,
    gene=None,
    annotation=None,
    coord="umap",
    vmin="p1",
    vmax="p99",
    cmap="jet",
    title=None,
    width=900,
    height=750,
    colors=None,
    palette=None,
    size=None,
    show=True,
    downsample=None,
    target_fill=0.05,
    normalize_per_cell=True,
    clip_norm_value=10,
    renderer="notebook",
):
    """
    Plot interactive embedding plot with plotly for a given AnnData object or path of .h5ad.

    Parameters
    ----------
    adata : _type_
            _description_
    obs : _type_, optional
            _description_, by default None
    variable : _type_, optional
            _description_, by default None
    gene : _type_, optional
            _description_, by default None
    annotation: str
            column to pass to mouse over
    coord : str, optional
            _description_, by default "umap"
    vmin : str, optional
            _description_, by default 'p1'
    vmax : str, optional
            _description_, by default 'p99'
    cmap : str, optional
            _description_, by default 'jet'
    title : _type_, optional
            _description_, by default None
    width : int, optional
            _description_, by default 1000
    height : int, optional
            _description_, by default 800
    colors : _type_, optional
            _description_, by default None
    palette : _type_, optional
            _description_, by default None
    size : _type_, optional
            _description_, by default None
    target_fill : float, optional
            _description_, by default 0.05
    show : bool, optional
            _description_, by default True
    renderer : str, optional
            Plotly renderer backend, by default ``"notebook"``.
            Available renderers: ``'plotly_mimetype'``, ``'jupyterlab'``,
            ``'nteract'``, ``'vscode'``, ``'notebook'``,
            ``'notebook_connected'``, ``'kaggle'``, ``'azure'``, ``'colab'``,
            ``'cocalc'``, ``'databricks'``, ``'json'``, ``'png'``, ``'jpeg'``,
            ``'jpg'``, ``'svg'``, ``'pdf'``, ``'browser'``, ``'firefox'``,
            ``'chrome'``, ``'chromium'``, ``'iframe'``,
            ``'iframe_connected'``, ``'sphinx_gallery'``,
            ``'sphinx_gallery_png'``.

    Returns
    -------
    _type_
            _description_
    """
    import plotly.express as px
    import plotly.io as pio

    if renderer is not None:
        pio.renderers.default = renderer
    use_col = gene if gene is not None else variable
    if gene is not None:
        assert adata is not None, "`gene` provided, `adata` must be provided too."
    if adata is not None:
        adata = load_adata(adata)
    use_adata = None
    if gene is not None:  # adata is not None
        if adata.isbacked:  # type: ignore
            use_adata = adata[:, gene].to_memory()  # type: ignore
        else:
            use_adata = adata[:, gene].copy()  # type: ignore
        if normalize_per_cell:
            use_adata = normalize_mc_by_cell(
                use_adata=use_adata,
                normalize_per_cell=normalize_per_cell,
                clip_norm_value=clip_norm_value,
                hypo_score=False,
            )
    else:
        if adata is not None:
            use_adata = adata  # backed
    if obs is None and use_adata is None:
        raise ValueError("Either `adata` or `obs` must be provided.")
    if obs is None:
        obs = use_adata.obs.copy()  # type: ignore
    else:  # obs not none
        obs = load_obs(obs)
        if use_adata is not None:
            overlap_idx = obs.index.intersection(use_adata.obs_names)
            obs = obs.loc[overlap_idx]
            use_adata = use_adata[overlap_idx, :]  # type: ignore

    if gene is not None:
        obs[gene] = use_adata.to_df()[gene].values  # type: ignore
    # cols=set(obs.columns.tolist())
    if f"X_{coord}" in use_adata.obsm:  # type: ignore
        # print(use_adata.obsm[f'X_{coord}'])
        obs[f"{coord}_0"] = use_adata.obsm[f"X_{coord}"][:, 0]  # type: ignore
        obs[f"{coord}_1"] = use_adata.obsm[f"X_{coord}"][:, 1]  # type: ignore
        # print(obs.head())
    if adata is not None and adata.isbacked:  # type: ignore
        adata.file.close()  # type: ignore
    # downsample obs for large dataset
    n_points = obs.shape[0]
    if downsample is not None and n_points > downsample:
        sample_idx = np.random.choice(
            n_points, size=downsample, replace=False
        )  # numbers
        obs = obs.iloc[sample_idx]

    if obs.dtypes[use_col] not in ["object", "category"]:  # gene
        vmin_quantile = float(int(vmin.replace("p", "")) / 100)
        vmax_quantile = float(int(vmax.replace("p", "")) / 100)
        # print(vmin_quantile,vmax_quantile,obs[use_col],obs.dtypes[use_col])
        range_color = [
            obs[use_col].quantile(vmin_quantile),
            obs[use_col].quantile(vmax_quantile),
        ]
        color_discrete_map = None
        category_orders = None
    else:  # categorical variable
        if colors is None:
            # color_discrete_map=get_colors(use_adata,use_col,palette=palette)
            color_discrete_map = load_color_palette(
                palette=palette, adata=use_adata, groups=use_col
            )
        else:
            color_discrete_map = colors
        if color_discrete_map is not None:
            keys = list(color_discrete_map.keys())  # type: ignore
            for k in keys:
                if k not in obs[use_col].unique().tolist():
                    del color_discrete_map[k]  # type: ignore
        range_color = None
        order = sorted(obs[use_col].dropna().unique().tolist())
        obs[use_col] = pd.Categorical(obs[use_col], categories=order, ordered=True)
        category_orders = {use_col: order}
    keep_cols = ["cell", f"{coord}_0", f"{coord}_1"]
    if variable is not None:
        keep_cols.append(variable)
    if gene is not None:
        keep_cols.append(gene)
    if annotation is not None:
        keep_cols.append(annotation)
    obs = obs.reset_index(names="cell").loc[:, keep_cols]

    # Create Plotly interactive scatter plot
    hover_data = {  # Fields to show on hover
        "cell": True,  # cell ID
        f"{coord}_0": ":0.3f",  # UMAP coordinates rounded to 3 decimals
        f"{coord}_1": ":0.3f",
    }
    if variable is not None:
        hover_data[variable] = True  # type: ignore # when plotting gene expression, also show cell types when mouse hover
    if annotation is not None:
        hover_data[annotation] = True  # type: ignore # when plotting gene expression, also show cell types when mouse hover
    if gene is not None:
        hover_data[gene] = ":.3f"  # type: ignore
    fig = px.scatter(
        obs,
        x=f"{coord}_0",  # UMAP first dimension → X axis
        y=f"{coord}_1",  # UMAP second dimension → Y axis
        color=use_col,
        category_orders=category_orders,
        hover_data=hover_data,
        range_color=range_color,
        color_discrete_sequence=px.colors.qualitative.D3,  # color palette (professional, unobtrusive)
        color_discrete_map=color_discrete_map,
        color_continuous_scale=cmap,  # ["blue", "red"],
        template="plotly_white",
        render_mode="webgl",  # use WebGL rendering for better performance with large datasets
    )
    fig.update_xaxes(
        range=[obs[f"{coord}_0"].min() - 0.5, obs[f"{coord}_0"].max() + 0.5],
        tickfont_size=12,
    )
    fig.update_yaxes(
        range=[obs[f"{coord}_1"].min() - 0.5, obs[f"{coord}_1"].max() + 0.5],
        tickfont_size=12,
    )

    if size is None:
        # Blend an area-based marker estimate with a log-based fallback so total point count and canvas size both matter.
        # Increased target_fill and scaling to make markers bigger
        marker_diam_area = 2 * np.sqrt(
            (width * height * target_fill) / (np.pi * n_points)
        )
        marker_diam_log = 16 - 2 * np.log10(n_points)
        marker_diam = 0.7 * marker_diam_area + 0.5 * marker_diam_log
        size = int(np.clip(marker_diam, 4, 20))
    if n_points < 500000:
        opacity = 0.8
    else:
        opacity = 0.6
    # logger.debug(f"{variable},{gene},{use_col}")
    # print(color_discrete_map,size,opacity)
    fig.update_traces(
        marker=dict(size=size, opacity=opacity, line=dict(width=0.12, color="black")),
        selector=dict(mode="markers"),
    )
    if title is None:
        title = f"{coord.upper()} Visualization (Colored by {use_col})"
    fig.update_layout(
        title=dict(
            text=title,
            font_size=16,
            x=0.5,  # center the title
            pad=dict(t=10),
        ),
        xaxis_title=f"{coord}_0".upper(),
        yaxis_title=f"{coord}_1".upper(),
        autosize=True,
        width=width,
        height=height,
        legend_title=use_col,  # legend title
        legend=dict(
            font_size=12,
            itemsizing="constant",  # important: fix legend marker size so it's not affected by scatter points
            itemwidth=30,
            borderwidth=0.1,  # legend item width; larger value increases the marker size
        ),
    )
    if show:
        filename = f"{coord}.{use_col}"
        show_fig(fig, filename=filename)
    else:
        return fig
    # html=fig2div(fig,filename='umap_plot')
    # return HttpResponse(html)


def show_fig(fig, filename="plot"):
    """Display a Plotly figure with an interactive toolbar configuration.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to display.
    filename : str, optional
        Base filename used when exporting the figure to SVG via the toolbar
        download button. Default is ``'plot'``.
    """
    interactive_config = {
        "displayModeBar": "hover",
        "showLink": False,
        "linkText": "Edit on plotly",
        "scrollZoom": True,
        "displaylogo": False,
        "toImageButtonOptions": {"format": "svg", "filename": filename},
        "modeBarButtonsToRemove": [
            "sendDataToCloud"
        ],  # 'zoomIn2d','zoomOut2d','zoom2d','zoom3d','pan2d'
        "editable": True,
        "autosizable": True,  #'responsive':True, 'fillFrame':True,
        "edits": {
            "titleText": True,
            "legendPosition": True,
            "colorbarTitleText": True,
            "shapePosition": True,
            "annotationPosition": True,
            "annotationText": True,
            "axisTitleText": True,
            "legendText": True,
            "colorbarPosition": True,
        },
    }
    fig.show(config=interactive_config)


def plot_categorical(
    adata,
    ax=None,
    basis="umap",
    groupby="Subclass",
    coding=True,
    coded_marker=True,
    save=None,
    palette=None,
    sheet_name=None,
    show=True,
    figsize=(4, 3.5),
    ncol=None,
    fontsize=5,
    legend_fontsize=5,
    legend_kws=None,
    legend_title_fontsize=5,
    marker_fontsize=4,
    marker_pad=0.1,
    linewidth=0.5,
    axis_format="tiny",
    alpha=0.7,
    text_kws=None,
    **kwargs,
):
    """Create a categorical scatter plot (e.g. UMAP) colored by a grouping variable.

    Parameters
    ----------
    adata : AnnData or str
            AnnData object or path to an ``.h5ad`` file.
    ax : matplotlib.axes.Axes or None, optional
            Pre-existing axes to draw on. If ``None``, a new figure is created.
    basis : str, optional
            Embedding key in ``adata.obsm`` (without the ``X_`` prefix).
            Default is ``'umap'``.
    groupby : str, optional
            Column in ``adata.obs`` used to color the scatter plot.
            Default is ``'Subclass'``.
    coding : bool, optional
            If ``True``, assign numeric codes to each group and display them as
            markers on the plot. Default is ``True``.
    coded_marker : bool, optional
            If ``True``, show the numeric code inside each group's marker.
            Default is ``True``.
    save : str or None, optional
            Path to save the figure. If ``None``, the figure is not saved.
    palette : str, dict, or None, optional
            Color palette. Can be a path to an Excel palette file, a
            ``{category: hex_color}`` dict, or ``None`` to use colours stored
            in ``adata.uns``.
    sheet_name : str or None, optional
            Sheet name to read from the palette Excel file. Defaults to
            *groupby*.
    show : bool, optional
            Whether to call ``plt.show()``. Default is ``True``.
    figsize : tuple of float, optional
            Figure size ``(width, height)`` in inches. Default is ``(4, 3.5)``.
    ncol : int or None, optional
            Number of columns in the legend. ``None`` lets matplotlib decide.
    fontsize : int, optional
            Font size for text annotations on the plot. Default is ``5``.
    legend_fontsize : int, optional
            Font size for legend entries. Default is ``5``.
    legend_kws : dict or None, optional
            Extra keyword arguments forwarded to the legend.
    legend_title_fontsize : int, optional
            Font size for the legend title. Default is ``5``.
    marker_fontsize : int, optional
            Font size for coded markers. Default is ``4``.
    marker_pad : float, optional
            Padding around coded markers. Default is ``0.1``.
    linewidth : float, optional
            Edge linewidth of scatter points. Default is ``0.5``.
    axis_format : str, optional
            Axis appearance style (e.g. ``'tiny'``, ``'arrow'``). Default is
            ``'tiny'``.
    alpha : float, optional
            Point transparency. Default is ``0.7``.
    text_kws : dict or None, optional
            Extra keyword arguments for text annotations.
    **kwargs
            Additional keyword arguments (such as outline, outline_kws et.al) passed to
            :func:`~adataviz.utils.categorical_scatter`.
    """
    from pandas.api.types import is_categorical_dtype

    if basis.startswith("X_"):
        basis = basis.replace("X_", "")
    if sheet_name is None:
        sheet_name = groupby
    adata = load_adata(adata)
    if not is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype("category")
    colors = None
    if palette is not None:
        if isinstance(palette, dict):
            colors = palette.copy()
        elif isinstance(palette, str):
            try:
                colors = pd.read_excel(
                    os.path.expanduser(palette), sheet_name=sheet_name, index_col=0
                ).Hex.to_dict()
                keys = list(colors.keys())
                existed_vals = adata.obs[groupby].unique().tolist()
                for k in existed_vals:
                    if k not in keys:
                        colors[k] = "gray"
                for k in keys:
                    if k not in existed_vals:
                        del colors[k]
            except Exception:
                colors = None
    if colors is None:
        if f"{groupby}_colors" in adata.uns:
            colors = {
                cluster: color
                for cluster, color in zip(
                    adata.obs[groupby].cat.categories.tolist(),
                    adata.uns[f"{groupby}_colors"],
                )
            }
        else:
            colors = _assign_default_colors(adata, groupby)
    else:
        adata.uns[groupby + "_colors"] = [
            colors.get(k, "grey") for k in adata.obs[groupby].cat.categories.tolist()
        ]

    hue = groupby
    text_anno = groupby
    text_kws = {} if text_kws is None else text_kws
    text_kws.setdefault("fontsize", fontsize)
    kwargs.setdefault("hue", hue)
    kwargs.setdefault("text_anno", text_anno)
    kwargs.setdefault("text_kws", text_kws)
    kwargs.setdefault("luminance", 0.65)
    kwargs.setdefault("dodge_text", False)
    kwargs.setdefault("axis_format", axis_format)
    kwargs.setdefault("show_legend", True)
    kwargs.setdefault("marker_fontsize", marker_fontsize)
    kwargs.setdefault("marker_pad", marker_pad)
    kwargs.setdefault("linewidth", linewidth)
    kwargs.setdefault("alpha", alpha)
    kwargs["coding"] = coding
    kwargs["coded_marker"] = coded_marker
    legend_kws = {} if legend_kws is None else legend_kws
    default_lgd_kws = dict(
        fontsize=legend_fontsize, title=groupby, title_fontsize=legend_title_fontsize
    )
    if ncol is not None:
        default_lgd_kws["ncol"] = ncol
    for k in default_lgd_kws:
        legend_kws.setdefault(k, default_lgd_kws[k])
    kwargs.setdefault(
        "dodge_kws",
        {
            "arrowprops": {
                "arrowstyle": "->",
                "fc": "grey",
                "ec": "none",
                "connectionstyle": "angle,angleA=-90,angleB=180,rad=5",
            },
            "autoalign": "xy",
        },
    )
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=300)
    categorical_scatter(
        data=adata[adata.obs[groupby].notna(),],
        ax=ax,
        basis=basis,
        palette=colors,
        legend_kws=legend_kws,
        **kwargs,
    )

    if save is not None and save is not False:
        plt.savefig(os.path.expanduser(save), bbox_inches="tight")  # ,,dpi=300
    if show:
        plt.show()


def plot_gene(
    adata,
    obs=None,
    groupby=None,
    gene="CADM1",
    query_str=None,
    title=None,
    palette=None,
    hue_norm=None,
    cbar_kws=dict(extendfrac=0.1),
    axis_format="tiny",
    scatter_kws={},
    obsm=None,
    basis="umap",
    normalize_per_cell=False,
    stripplot=False,
    hypo_score=False,
    ylim=None,
    clip_norm_value=10,
    min_cells=3,
    cmap="parula",
    prefix=None,
):
    """Plot gene expression on an embedding and optionally a violin/strip plot.

    Reads the gene from a backed ``.h5ad`` file, projects expression values
    onto the chosen embedding coordinates, and saves the scatter plot to PDF.
    When *groupby* is provided, a violin plot grouped by that column is also
    generated.

    Parameters
    ----------
    adata : str
            Path to an ``.h5ad`` file.
    obs : str, DataFrame, or None, optional
            Cell metadata. Can be a path to a TSV file or a DataFrame.
            If ``None``, uses ``adata.obs``.
    groupby : str or None, optional
            Column in *obs* to group cells for the violin/strip plot.
    gene : str, optional
            Gene name to visualise. Default is ``'CADM1'``.
    query_str : str or None, optional
            Pandas query string to subset *obs* before plotting.
    title : str or None, optional
            Plot title. Defaults to *query_str*, *groupby*, or *gene*.
    palette : str, dict, or None, optional
            Color palette for the violin plot. Can be a path to an Excel
            palette file, a ``{group: hex_color}`` dict, or a column name in
            ``adata.obs`` containing hex colours.
    hue_norm : tuple or None, optional
            ``(vmin, vmax)`` to normalise the colour mapping for the scatter.
    cbar_kws : dict, optional
            Keyword arguments for the colourbar. Default is
            ``dict(extendfrac=0.1)``.
    axis_format : str, optional
            Axis appearance style. Default is ``'tiny'``.
    scatter_kws : dict, optional
            Extra keyword arguments forwarded to
            :func:`~adataviz.utils.continuous_scatter`.
    obsm : str, AnnData, or None, optional
            AnnData (or path) supplying alternative ``obsm`` embeddings and
            extra ``obs`` columns.
    basis : str, optional
            Embedding key (without ``X_`` prefix). Default is ``'umap'``.
    normalize_per_cell : bool, optional
            Whether to normalise methylation values per cell. Default is
            ``False``.
    stripplot : bool, optional
            If ``True``, overlay a strip plot on the violin. Default is
            ``False``.
    hypo_score : bool, optional
            Whether to compute hypomethylation score. Default is ``False``.
    ylim : tuple or None, optional
            Y-axis limits for the violin plot.
    clip_norm_value : float, optional
            Clipping value for per-cell normalisation. Default is ``10``.
    min_cells : int, optional
            Minimum number of cells required for a group to be included in
            the violin plot. Default is ``3``.
    cmap : str, optional
            Colourmap for the scatter plot. Default is ``'parula'``.
    prefix : str or None, optional
            File-name prefix for saved PDFs. Defaults to
            ``'{title}.{gene}.{groupby}'``.
    """
    import seaborn as sns

    # sc.set_figure_params(dpi=100,dpi_save=300,frameon=False)
    if title is None:
        if query_str is not None:
            title = query_str
        else:
            title = groupby if groupby is not None else gene
    raw_adata = anndata.read_h5ad(os.path.expanduser(adata), backed="r")
    adata = raw_adata[:, gene].to_memory()
    raw_adata.file.close()  # close the file to save memory
    if normalize_per_cell:
        adata = normalize_mc_by_cell(
            use_adata=adata,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
            hypo_score=hypo_score,
        )
    is_open = False
    if obsm is not None:
        if isinstance(obsm, str):
            obsm = anndata.read_h5ad(os.path.expanduser(obsm), backed="r")
            is_open = True
        assert isinstance(obsm, anndata.AnnData), (
            "obsm should be an anndata object or a path to an h5ad file."
        )
        keep_cells = list(set(adata.obs_names.tolist()) & set(obsm.obs_names.tolist()))
        adata = adata[keep_cells, :]
        adata.obsm = obsm[keep_cells].obsm
        cur_cols = adata.obs.columns.tolist()
        for col in obsm.obs.columns.tolist():
            if col not in cur_cols:
                adata.obs[col] = obsm.obs.loc[adata.obs_names, col].tolist()
    if is_open:
        obsm.file.close()
    if obs is not None:
        if isinstance(obs, str):
            obs = pd.read_csv(os.path.expanduser(obs), sep="\t", index_col=0)
        else:
            obs = obs.copy()
    else:
        obs = adata.obs.copy()
    if query_str is not None:
        obs = obs.query(query_str)
    overlapped_cells = list(set(adata.obs_names.tolist()) & set(obs.index.tolist()))
    obs = obs.loc[overlapped_cells]
    adata = adata[overlapped_cells, :]  # type: ignore
    adata.obs = obs.loc[adata.obs_names.tolist()]
    # print(adata.shape)
    # read color palette
    if groupby is not None and palette is not None:
        if isinstance(palette, dict):
            color_palette = palette.copy()
        elif isinstance(palette, str) and os.path.exists(os.path.expanduser(palette)):
            palette = os.path.abspath(os.path.expanduser(palette))
            D = pd.read_excel(palette, sheet_name=None, index_col=0)
            color_palette = D[groupby].Hex.to_dict()
        elif isinstance(palette, str):
            color_palette = (
                adata.obs.reset_index()
                .loc[:, [groupby, palette]]
                .drop_duplicates()
                .dropna()
                .set_index(groupby)[palette]
                .to_dict()
            )
        else:
            color_palette = None
    else:
        color_palette = None
    # plot gene on given cordinate base
    # fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    # output=os.path.join(figdir, f"{title}.{gene}.{basis}.pdf")
    # sc.pl.embedding(adata, basis=basis,
    # 				wspace=0.1, color=[gene],use_raw=False,
    # 				ncols=2, vmin='p5', vmax='p95', frameon=False,
    # 				show=False,cmap=cmap,ax=ax)
    # colorbar = fig.axes[-1]
    # cur_pos=colorbar.get_position()
    # colorbar.set_position([cur_pos.x0,(1-cur_pos.height/2)/2,cur_pos.width, cur_pos.height / 2])
    # fig.savefig(output) # transparent=True,bbox_inches='tight',dpi=300

    if prefix is None:
        prefix = f"{title}.{gene}.{groupby}"
    adata.obs[gene] = adata.to_df().loc[adata.obs_names.tolist(), gene].tolist()
    # print(hue_norm)
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    continuous_scatter(
        data=adata,
        ax=ax,
        cmap=cmap,
        hue_norm=hue_norm,
        cbar_kws=cbar_kws,
        hue=gene,
        axis_format=axis_format,
        text_anno=None,
        basis=basis,
        **scatter_kws,
    )
    fig.savefig(
        f"{prefix}.{basis}.pdf", bbox_inches="tight"
    )  # transparent=True,bbox_inches='tight',dpi=300

    if groupby is not None:
        # boxplot
        data = adata.to_df()
        data[groupby] = adata.obs.loc[data.index.tolist(), groupby].tolist()
        vc = data[groupby].value_counts()
        N = vc.shape[0]
        if color_palette is not None:
            keep_groups = list(
                set(list(color_palette.keys()))
                & set(vc[vc >= min_cells].index.tolist())
            )
            data = data.loc[data[groupby].isin(keep_groups)]
        vc = vc.to_dict()
        order = data.groupby(groupby)[gene].median().sort_values().index.tolist()
        width = max(5, N * 0.5)
        plt.figure(figsize=(width, 3.5))
        if stripplot:
            ax = sns.stripplot(
                data=data,
                jitter=0.4,
                edgecolor="white",
                x=groupby,
                y=gene,
                palette=color_palette,
                order=order,
                size=0.5,
            )
        else:
            ax = None
        # ax = sns.boxplot(data=data, x=groupby, y=gene, palette=color_palette, ax=ax,  # hue=groupby,
        # 				fliersize=0.5, notch=False, showfliers=False, saturation=0.6, order=order)
        # boxplot are incorrect for some cases when there are many 0, median and lower quartile are often at zero; use violinplot in stead.
        ax = sns.violinplot(
            data=data,
            x=groupby,
            y=gene,
            palette=color_palette,
            ax=ax,  # hue=groupby,
            saturation=0.6,
            order=order,
            density_norm="width",
            cut=0,
            bw_adjust=0.5,
        )
        # ax=sns.swarmplot(data=data,palette=color_palette,\
        #                   edgecolor='white',x=groupby,y=gene,\
        #                   order=order)
        if ylim is not None:
            ax.set_ybound(ylim)
        ax.set_xticklabels([f"{label} ({vc[label]})" for label in order])
        title = title.replace(" ", ".")
        ax.set_title(title)
        ax.xaxis.label.set_visible(False)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=-45, ha="left")
        plt.savefig(f"{prefix}.boxplot.pdf", bbox_inches="tight")


def stacked_barplot(
    obs,
    groupby="Age",
    column="CellClass",
    x_order=None,
    y_order=None,
    linewidth=0.1,
    palette=None,
    width=None,
    height=None,
    xticklabels_kws=None,
    save=False,
    lgd_kws=None,
    gap=0.05,
    sort_by=None,
):
    """Plot stacked barplot to show cell type composition in each *groupby* group.

    For example Age, brain regions, etc.

    Examples::

            stacked_barplot(column='MajorType', width=3.5, height=6)
            stacked_barplot(column='CellClass', width=3.5, height=3)
    """
    if isinstance(obs, pd.DataFrame):
        data = obs.copy()
    elif isinstance(obs, str) and obs.endswith(".h5ad"):
        obs_path = os.path.abspath(os.path.expanduser(obs))
        adata = anndata.read_h5ad(obs_path, backed="r")
        data = adata.obs
        del adata
    elif obs.endswith(".csv"):
        obs_path = os.path.abspath(os.path.expanduser(obs))
        data = pd.read_csv(obs_path, index_col=0)
    else:
        obs_path = os.path.abspath(os.path.expanduser(obs))
        data = pd.read_csv(obs_path, index_col=0, sep="\t")
    xticklabels_kws = {} if xticklabels_kws is None else xticklabels_kws
    xticklabels_kws.setdefault("rotation", -45)
    xticklabels_kws.setdefault("rotation_mode", "anchor")
    xticklabels_kws.setdefault(
        "horizontalalignment", "left"
    )  # see ?matplotlib.axes.Axes.tick_params
    xticklabels_kws.setdefault("verticalalignment", "center")
    if palette is not None:
        if isinstance(palette, dict):
            color_palette = palette.copy()
        elif isinstance(palette, str) and os.path.exists(os.path.expanduser(palette)):
            palette = os.path.abspath(os.path.expanduser(palette))
            D = pd.read_excel(palette, sheet_name=None, index_col=0)
            color_palette = D[column].Hex.to_dict()
            keys = list(color_palette.keys())
            for k in data[column].unique():
                if k not in keys:
                    color_palette[k] = "gray"
        else:
            color_palette = palette
    else:
        color_palette = None
    df = (
        data.groupby(groupby)[column].value_counts(normalize=True).unstack(level=column)
    )
    if sort_by is not None:
        df.sort_values(sort_by, ascending=True, inplace=True)
    else:
        if x_order is None:
            x_order = sorted(df.index.tolist())
        if y_order is None:
            y_order = sorted(df.columns.tolist())
        df = df.loc[x_order, y_order]
    if width is None:
        width = max(df.shape[0] * 0.45, 10)
        if width < 2.5:
            width = 2.5
    if height is None:
        height = max(df.shape[1] * 0.5, 8)
        if height < 3.5:
            height = 3.5
    plt.figure()
    ax = df.plot.bar(
        stacked=True,
        align="edge",
        width=1 - gap,
        edgecolor="black",
        linewidth=linewidth,
        figsize=(width, height),
        color=color_palette,
    )
    ax.set_xlim(0, df.shape[0])
    ax.set_ylim(0, 1)
    labels = [tick.get_text() for tick in ax.get_xticklabels()]
    ax.set_xticks(
        ticks=np.arange(0.5, df.shape[0], 1), labels=labels, **xticklabels_kws
    )  # ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0.5,df.shape[1],1))) #ticker.MultipleLocator(0.5)
    ax.xaxis.label.set_visible(False)
    ax.tick_params(
        axis="y",  # both
        which="both",
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
        top=False,
        labeltop=False,  # bottom=False,labelbottom=False
    )
    # ax.xaxis.set_tick_params(axis='x')
    lgd_kws = lgd_kws if lgd_kws is not None else {}  # bbox_to_anchor=(x,-0.05)
    lgd_kws.setdefault("frameon", True)
    lgd_kws.setdefault("ncol", 1)
    lgd_kws["loc"] = "upper left"
    lgd_kws.setdefault("borderpad", 0.1 * (1 / 25.4) * 72)  # 0.1mm
    lgd_kws.setdefault("markerscale", 1)
    lgd_kws.setdefault("handleheight", 1)  # font size, units is points
    lgd_kws.setdefault("handlelength", 1)  # font size, units is points
    lgd_kws.setdefault(
        "borderaxespad", 0.1
    )  # The pad between the axes and legend border, in font-size units.
    lgd_kws.setdefault(
        "handletextpad", 0.3
    )  # The pad between the legend handle and text, in font-size units.
    lgd_kws.setdefault(
        "labelspacing", 0.1
    )  # gap height between two Patches,  0.05*mm2inch*72
    lgd_kws.setdefault("columnspacing", 1)
    lgd_kws.setdefault("bbox_to_anchor", (1, 1))
    lgd_kws.setdefault("title", column)
    ax.legend(**lgd_kws)
    ax.grid(False)
    if save:
        outdir = os.path.dirname(os.path.expanduser(save))
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        plt.savefig(
            save, bbox_inches="tight"
        )  # transparent=True,bbox_inches='tight',dpi=300
    else:
        plt.show()
    return ax


def pieplot(obs, groupby="Age", palette=None, order=None, save=None, explode=0.05):
    """Create a pie chart showing the distribution of cells across groups.

    Parameters
    ----------
    obs : DataFrame, str
            Cell metadata. Can be a DataFrame, a path to an ``.h5ad`` file, a
            ``.csv`` file, or a tab-separated text file.
    groupby : str, optional
            Column in *obs* whose value counts define the pie slices.
            Default is ``'Age'``.
    palette : str, dict, or None, optional
            Color palette. Can be a ``{group: hex_color}`` dict or a path to
            an Excel palette file. ``None`` uses matplotlib defaults.
    order : list or None, optional
            Explicit ordering of pie slices. If ``None``, categories are
            sorted alphabetically.
    save : str or None, optional
            Path to save the figure. If ``None``, saves to
            ``'{groupby}.piechart.pdf'``.
    explode : float, optional
            Fraction by which each slice is offset from the centre.
            Default is ``0.05``.
    """
    # colors=None
    if isinstance(obs, pd.DataFrame):
        data = obs.copy()
    elif isinstance(obs, str) and obs.endswith(".h5ad"):
        obs_path = os.path.abspath(os.path.expanduser(obs))
        print(f"Reading adata: {obs}")
        adata = anndata.read_h5ad(obs_path, backed="r")
        # if f'{groupby}_colors' in adata.uns:
        #     colors={k:v for k,v in zip(adata.obs[groupby].cat.categories.tolist(),
        #                                adata.uns[f'{groupby}_colors'])}
        # else:
        #     colors=None
        data = adata.obs
        del adata
    elif obs.endswith(".csv"):
        obs_path = os.path.abspath(os.path.expanduser(obs))
        data = pd.read_csv(obs_path, index_col=0)
    else:
        obs_path = os.path.abspath(os.path.expanduser(obs))
        data = pd.read_csv(obs_path, index_col=0, sep="\t")

    if palette is not None:
        if isinstance(palette, dict):
            color_palette = palette.copy()
        else:
            palette = os.path.abspath(os.path.expanduser(palette))
            D = pd.read_excel(palette, sheet_name=None, index_col=0)
            color_palette = D[groupby].Hex.to_dict()
    else:
        color_palette = None
    D = data[groupby].value_counts()
    if order is None:
        order = list(sorted(D.keys()))

    plt.figure()
    plt.pie(
        [D[k] for k in order],
        labels=order,
        colors=[color_palette[k] for k in order] if color_palette else None,
        explode=[explode] * len(order),
        autopct="%.1f%%",
    )
    # Add title to the chart
    plt.title("Distribution of #of cells across different stages")
    if save is not None:
        output = os.path.abspath(os.path.expanduser(save))
    else:
        output = f"{groupby}.piechart.pdf"
    plt.savefig(
        output, bbox_inches="tight"
    )  # transparent=True,bbox_inches='tight',dpi=300
    plt.show()


def plot_pseudotime(
    pseudotime="dpt_pseudotime.tsv",
    groupby="Age",
    y="dpt_pseudotime",
    hue=None,
    figsize=(5, 3.5),
    save=None,
    rotate=None,
    ylabel="Pseudotime",
    palette=None,
):
    """Plot pseudotime as a violin plot grouped by a categorical variable.

    Examples::

            plot_pseudotime(figsize=(6, 3.5), groupby='MajorType', rotate=-45)
            plot_pseudotime(figsize=(3.5, 3), groupby='CellClass')
            plot_pseudotime(figsize=(3.5, 3), groupby='Age')

    Parameters
    ----------
    pseudotime : str, optional
            Path to a TSV file with pseudotime values. Default is
            ``'dpt_pseudotime.tsv'``.
    groupby : str, optional
            Column to group cells by. Default is ``'Age'``.
    y : str, optional
            Column containing pseudotime values. Default is
            ``'dpt_pseudotime'``.
    hue : str or None, optional
            Additional grouping variable for the violin plot.
    figsize : tuple of float, optional
            Figure size ``(width, height)``. Default is ``(5, 3.5)``.
    save : str or None, optional
            Path to save the figure. If ``None``, saves to
            ``'{groupby}.pseudotime_violin.pdf'``.
    rotate : float or None, optional
            Rotation angle for x-tick labels.
    yield : str, optional
            Y-axis label. Default is ``'Pseudotime'``.
    palette : str, dict, or None, optional
            Color palette. Can be a ``{group: hex_color}`` dict or a path
            to an Excel palette file.
    """
    import seaborn as sns

    if palette is not None:
        if isinstance(palette, dict):
            color_palette = palette.copy()
        else:
            palette = os.path.abspath(os.path.expanduser(palette))
            D = pd.read_excel(palette, sheet_name=None, index_col=0)
            color_palette = D[groupby].Hex.to_dict()
    else:
        color_palette = None

    data = pd.read_csv(os.path.expanduser(pseudotime), sep="\t", index_col=0)
    data.dpt_pseudotime.replace({np.inf: 1}, inplace=True)
    order = data.groupby(groupby)[y].mean().sort_values(ascending=True).index.tolist()
    if hue is not None:
        hue_order = (
            data.groupby(hue)[y].mean().sort_values(ascending=True).index.tolist()
        )
    else:
        hue_order = None
    plt.figure(figsize=figsize)
    # ax = sns.swarmplot(data=data, palette=color_palette, \
    # 				   edgecolor='white', x=groupby, y=y, \
    # 				   order=order)
    ax = sns.violinplot(
        data=data,
        x=groupby,
        y=y,
        scale="count",
        bw=0.2,
        inner=None,
        saturation=0.6,
        palette=color_palette,
        order=order,
        hue=hue,
        hue_order=hue_order,
    )
    # plt.legend(frameon=True)
    ax.set_ylabel(ylabel)
    if rotate is not None:
        plt.setp(
            ax.xaxis.get_majorticklabels(),
            rotation=rotate,
            rotation_mode="anchor",
            horizontalalignment="left",
        )
    if save is None:
        outname = (
            groupby + ".pseudotime_violin.pdf"
            if hue is None
            else groupby + f"_{hue}.pseudotime_violin.pdf"
        )
    else:
        outname = os.path.abspath(os.path.expanduser(save))
    plt.savefig(
        outname, bbox_inches="tight"
    )  # transparent=True,bbox_inches='tight',dpi=300
    plt.show()


def stacked_violinplot(
    adata,
    use_genes=None,
    groupby="Age",
    cell_groups=None,
    parent=None,
    figsize=(6, 4),
    cmap="viridis",
):
    """Plot a stacked violin plot for selected genes using scanpy.

    Parameters
    ----------
    adata : AnnData
            Annotated data matrix.
    use_genes : list of str or None, optional
            Gene names to plot.
    groupby : str, optional
            Column in ``adata.obs`` to group cells by. Default is ``'Age'``.
    cell_groups : list of str or None, optional
            List where the first element is an ``obs`` column used to subset
            cells (via *parent*).
    parent : str or None, optional
            Value within ``cell_groups[0]`` to select cells from.
    figsize : tuple of float, optional
            Figure size ``(width, height)``. Default is ``(6, 4)``.
    cmap : str, optional
            Colourmap. Default is ``'viridis'``.
    """
    import scanpy as sc

    ax = sc.pl.stacked_violin(
        adata[adata.obs[cell_groups[0]] == parent],
        var_names=use_genes,
        title=parent,
        colorbar_title="Avg mc frac",
        groupby=groupby,
        dendrogram=True,
        swap_axes=False,
        cmap=cmap,
        figsize=figsize,
        scale="count",
        standard_scale="obs",
        inner="quart",
        # stripplot=False,jitter=False,
        show=False,
        layer=None,
    )
    import matplotlib.ticker as ticker

    ax1 = ax["mainplot_ax"]
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax1.yaxis.set_tick_params(which="minor", left=True)
    ax1.grid(axis="y", linestyle="--", color="black", alpha=1, zorder=-5, which="minor")
    # plt.savefig(f"{fig_basename}.{groupby}.stacked_violin.pdf")
    plt.show()


def plot_genes(
    adata,
    query_str=None,
    obs=None,
    groupby="Subclass",
    parent_col=None,
    modality="RNA",  # mc or RNA
    use_raw=True,  # True for RNA
    expression_cutoff="p5",  # for RNA, could be int, median, mean of p5, p95 and so on
    genes=None,
    cell_type_order=None,
    gene_order=None,
    row_cluster=False,
    col_cluster=False,
    cmap="Greens_r",
    group_legend=False,
    parent_legend=False,
    title=None,
    palette=None,
    legend_kws=dict(extendfrac=0.1, extend="both", label="Mean mCG"),
    normalize_per_cell=True,
    clip_norm_value=10,
    hypo_score=False,
    figsize=(10, 3.5),
    marker="o",
    plot_kws={},
    transpose=False,
    outname="test.pdf",
):
    """Create a dot-heatmap for multiple genes using PyComplexHeatmap.

    Dot size encodes the fraction of cells expressing each gene (or showing
    hypomethylation), and colour encodes the mean expression (or
    methylation) level.

    Parameters
    ----------
    adata : str
            Path to an ``.h5ad`` file.
    query_str : str or None, optional
            Pandas query string to subset observations.
    obs : str, DataFrame, or None, optional
            Cell metadata. Path to a TSV/CSV file or a DataFrame. If ``None``,
            uses ``adata.obs``.
    groupby : str or list of str, optional
            Column(s) in *obs* to group cells. If a list, groups are joined
            with ``'+'``. Default is ``'Subclass'``.
    parent_col : str or None, optional
            Higher-level grouping column for hierarchical annotation.
    modality : str, optional
            Data modality (``'RNA'``, ``'ATAC'``, or methylation).
            Default is ``'RNA'``.
    use_raw : bool, optional
            Whether to use ``adata.raw`` for RNA. Default is ``True``.
    expression_cutoff : str or float, optional
            Threshold to determine whether a gene is "expressed". Can be a
            percentile string (``'p5'``), ``'median'``, ``'mean'``, or a
            numeric value. Default is ``'p5'``.
    genes : list of str
            Gene names to include.
    cell_type_order : list of str or None, optional
            Explicit ordering of cell types along the cell-type axis.
    gene_order : list of str or None, optional
            Explicit ordering of genes along the gene axis.
    row_cluster : bool, optional
            Whether to cluster rows. Default is ``False``.
    col_cluster : bool, optional
            Whether to cluster columns. Default is ``False``.
    cmap : str, optional
            Colourmap for the dot fill. Default is ``'Greens_r'``.
    group_legend : bool, optional
            Whether to show a legend for the groupby annotation.
            Default is ``False``.
    parent_legend : bool, optional
            Whether to show a legend for the parent annotation.
            Default is ``False``.
    title : str or None, optional
            Plot title. Defaults to *query_str* or *groupby*.
    palette : str, dict, or None, optional
            Color palette. Can be a path to an Excel file, a ``{group: hex}``
            dict, or a column name in ``adata.obs`` with hex colours.
    legend_kws : dict, optional
            Keyword arguments for the colourbar legend.
    normalize_per_cell : bool, optional
            Whether to normalise methylation per cell. Default is ``True``.
    clip_norm_value : float, optional
            Clipping value for per-cell normalisation. Default is ``10``.
    hypo_score : bool, optional
            Whether to compute hypomethylation score. Default is ``False``.
    figsize : tuple of float, optional
            Figure size ``(width, height)``. Default is ``(10, 3.5)``.
    marker : str, optional
            Marker shape for dots. Default is ``'o'``.
    plot_kws : dict, optional
            Extra keyword arguments forwarded to
            ``DotClustermapPlotter``.
    transpose : bool, optional
            If ``True``, genes are on columns and cell types on rows.
            Default is ``False``.
    outname : str, optional
            Output PDF file name. Default is ``'test.pdf'``.

    Returns
    -------
    tuple
            ``(plot_data, df_cols, cm1)`` – the tidy plotting DataFrame,
            the column-annotation DataFrame, and the
            ``DotClustermapPlotter`` object.
    """
    from PyComplexHeatmap import (
        HeatmapAnnotation,
        anno_label,
        anno_simple,
        DotClustermapPlotter,
    )

    assert genes is not None, "Please provide genes to plot."
    # adata could be single cell level or pseudobulk level (adata.layers['frac'] should be existed)
    raw_adata = anndata.read_h5ad(os.path.expanduser(adata), backed="r")

    all_vars = set(raw_adata.var_names.tolist())
    keep_genes = list(
        set(all_vars) & set(genes)
    )  # keep_genes=[g for g in all_vars if g in genes]
    error_genes = [g for g in genes if g not in keep_genes]
    if len(error_genes) > 0:
        print(f"genes not found in adata: {error_genes}")
    adata = raw_adata[:, keep_genes].to_memory()  # type: ignore
    if use_raw and adata.raw is not None:
        adata_raw = adata.raw[:, adata.var_names.tolist()].to_adata()
        adata.X = adata_raw[adata.obs_names.tolist(), adata.var_names.tolist()].X.copy()  # type: ignore
        del adata_raw
    raw_adata.file.close()  # close the file to save memory

    if obs is not None:
        if isinstance(obs, str):
            obs = pd.read_csv(os.path.expanduser(obs), sep="\t", index_col=0)
        else:
            obs = obs.copy()
    else:
        obs = adata.obs.copy()
    if query_str is not None:
        obs = obs.query(query_str)
    overlapped_cells = list(set(adata.obs_names.tolist()) & set(obs.index.tolist()))
    obs = obs.loc[overlapped_cells]
    adata = adata[overlapped_cells, :]  # type: ignore
    if isinstance(groupby, list):
        groupby1 = "+".join(groupby)
        obs[groupby1] = obs.loc[:, groupby].apply(
            lambda x: "+".join(x.astype(str).tolist()), axis=1
        )
        groupby = groupby1
    adata.obs[groupby] = obs.loc[adata.obs_names.tolist(), groupby].tolist()
    if title is None:
        if query_str is not None:
            title = query_str
        else:
            title = groupby if groupby is not None else "-".join(genes)
    if parent_col is not None and parent_col not in adata.obs.columns.tolist():
        adata.obs[parent_col] = obs.loc[adata.obs_names.tolist(), parent_col].tolist()

    if modality not in ["RNA", "ATAC"] and normalize_per_cell:
        adata = normalize_mc_by_cell(
            use_adata=adata,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
            hypo_score=hypo_score,
        )
    print(adata.shape)

    # read color palette
    color_palette = {}
    if palette is not None:
        if isinstance(palette, dict):
            color_palette[groupby] = palette.copy()
        elif isinstance(palette, str) and os.path.exists(os.path.expanduser(palette)):
            palette = os.path.abspath(os.path.expanduser(palette))
            D = pd.read_excel(palette, sheet_name=None, index_col=0)
            if groupby in D:
                color_palette[groupby] = D[groupby].Hex.to_dict()
            else:
                assert "+" in groupby, f"{groupby} not found in the palette file."
                for group in groupby.split("+"):
                    assert group in D, f"{group} not found in the palette file."
                    color_palette[group] = D[group].Hex.to_dict()
            if parent_col is not None:
                color_palette[parent_col] = D[parent_col].Hex.to_dict()
        elif isinstance(palette, str):
            color_palette[groupby] = (
                adata.obs.reset_index()
                .loc[:, [groupby, palette]]
                .drop_duplicates()
                .dropna()
                .set_index(groupby)[palette]
                .to_dict()
            )
            if parent_col is not None:
                color_palette[parent_col] = (
                    adata.obs.reset_index()
                    .loc[:, [parent_col, palette]]
                    .drop_duplicates()
                    .dropna()
                    .set_index(parent_col)[palette]
                    .to_dict()
                )
    else:
        color_palette = None

    data = adata.to_df()  # rows are cells or cell types, columns are genes
    if modality in ["RNA", "ATAC"] and isinstance(expression_cutoff, str):
        if expression_cutoff == "median":
            cutoff = data.stack().median()
        elif expression_cutoff == "mean":
            cutoff = data.stack().mean()
        else:  # quantile, such as p5,p95
            f = float(expression_cutoff.replace("p", ""))
            cutoff = data.stack().quantile(f / 100)
        expression_cutoff = cutoff

    data[groupby] = adata.obs.loc[data.index.tolist(), groupby].tolist()
    if parent_col is not None and parent_col in adata.obs.columns.tolist():
        group2parent = (
            adata.obs.loc[:, [groupby, parent_col]]
            .drop_duplicates()
            .set_index(groupby)[parent_col]
            .to_dict()
        )
    plot_data = data.groupby(groupby).mean().stack().reset_index()
    plot_data.columns = [groupby, "Gene", "Mean"]
    if "frac" in adata.layers:
        D = adata.to_df(layer="frac").stack().to_dict()
    else:
        if modality not in ["RNA", "ATAC"]:  # methylation, cutoff = 1
            assert normalize_per_cell, "Normalized methylation fraction is required"
            hypo_frac = data.groupby(groupby).agg(
                lambda x: x[x < 1].shape[0] / x.shape[0]
            )  # fraction of cells showing hypomethylation for the corresponding genes
            D = hypo_frac.stack().to_dict()
        else:  # for RNA
            print(f"Using expression cutoff: {expression_cutoff}")
            frac = data.groupby(groupby).agg(
                lambda x: x[x > expression_cutoff].shape[0] / x.shape[0]
            )  # raw count > expression_cutoff means the gene is expressed
            D = frac.stack().to_dict()
    plot_data["frac"] = (
        plot_data.loc[:, [groupby, "Gene"]]
        .apply(lambda x: tuple(x.tolist()), axis=1)
        .map(D)
    )
    # plot_data

    df_cols = pd.DataFrame(
        list(sorted(adata.obs[groupby].unique().tolist())), columns=[groupby]
    )
    if parent_col is not None:
        df_cols[parent_col] = df_cols[groupby].map(group2parent)
        df_cols.sort_values([parent_col, groupby], inplace=True)
    df_cols.index = df_cols[groupby].tolist()
    if cell_type_order is not None:
        rows = [ct for ct in cell_type_order if ct in df_cols.index.tolist()]
        df_cols = df_cols.loc[rows]
    col_ha_dict = {}
    if "+" in groupby:
        individual_groups = groupby.split("+")
        for ig in individual_groups:
            df_cols[ig] = df_cols[groupby].apply(
                lambda x: x.split("+")[individual_groups.index(ig)]
            )
            group_colors = {}
            for k in df_cols[ig].unique().tolist():
                group_colors[k] = color_palette[ig][k]
            col_ha_dict[ig] = anno_simple(
                df_cols[ig],
                colors=group_colors,
                add_text=False,
                legend=group_legend,
                height=3,
                label=ig,
            )
    df_cols.dropna(inplace=True)
    # df_cols.head()
    if parent_col is not None:
        parent_colors = {}
        axis = (
            1 if not transpose else 0
        )  # 1 for vertical (col annotation), 0 for horizontal
        for k in df_cols[parent_col].unique().tolist():
            parent_colors[k] = color_palette[parent_col][k]
        if "+" not in groupby:
            group_colors = {}
            for k in df_cols[groupby].unique().tolist():
                group_colors[k] = color_palette[groupby][k]
            col_ha = HeatmapAnnotation(
                axis=axis,
                label=anno_label(
                    df_cols[groupby],
                    colors=group_colors,
                    merge=True,
                    rotation=45,
                    fontsize=12,
                    arrowprops=dict(visible=False),
                ),
                group=anno_simple(
                    df_cols[groupby],
                    colors=group_colors,
                    add_text=False,
                    legend=group_legend,
                    height=3,
                    label=groupby,
                ),
                parent=anno_simple(
                    df_cols[parent_col],
                    colors=parent_colors,
                    add_text=False,
                    legend=parent_legend,
                    height=3,
                    label=parent_col,
                ),
            )
        else:
            col_ha_dict[parent_col] = anno_simple(
                df_cols[parent_col],
                colors=parent_colors,
                add_text=False,
                legend=parent_legend,
                height=3,
                label=parent_col,
            )
            col_ha = HeatmapAnnotation(**col_ha_dict, axis=axis, verbose=0)
        colnames = False

    else:
        axis = (
            1 if not transpose else 0
        )  # 1 for vertical (col annotation), 0 for horizontal
        if "+" not in groupby:
            group_colors = {}
            for k in df_cols[groupby].unique().tolist():
                group_colors[k] = color_palette[groupby][k]
            col_ha = HeatmapAnnotation(
                axis=axis,
                group=anno_simple(
                    df_cols[groupby],
                    colors=group_colors,
                    add_text=False,
                    legend=group_legend,
                    height=3,
                    label=groupby,
                ),
            )
        else:
            col_ha = HeatmapAnnotation(**col_ha_dict, axis=axis, verbose=0)
        colnames = True
    if not transpose:
        top_annotation = col_ha
        left_annotation = None
        x = groupby
        y = "Gene"
        x_order = df_cols.index.tolist()
        y_order = gene_order
        show_colnames = colnames
        show_rownames = True
    else:
        top_annotation = None
        left_annotation = col_ha
        y = groupby
        x = "Gene"
        y_order = df_cols.index.tolist()
        x_order = gene_order
        show_rownames = colnames
        show_colnames = True

    default_plot_kws = dict(
        marker=marker,
        grid=None,
        legend_vgap=8,
        dot_legend_marker=marker,
        cmap_legend_kws=legend_kws,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        row_cluster_method="ward",
        row_cluster_metric="euclidean",
        col_cluster_method="ward",
        col_cluster_metric="euclidean",
        col_names_side="top",
        row_names_side="left",
        show_rownames=show_rownames,
        show_colnames=show_colnames,
        row_dendrogram=False,
        # vmin=0,vmax=1.5,
        xticklabels_kws={
            "labelrotation": 45,
            "labelcolor": "blue",
            "labelsize": 10,
            "top": True,
        },
        yticklabels_kws={"labelcolor": "blue", "labelsize": 10, "left": True},
        spines=False,
    )
    for k in default_plot_kws:
        if k not in plot_kws:
            plot_kws[k] = default_plot_kws[k]

    plt.figure(figsize=figsize)
    cm1 = DotClustermapPlotter(
        data=plot_data,
        top_annotation=top_annotation,
        left_annotation=left_annotation,
        x_order=x_order,
        y_order=y_order,
        x=x,
        y=y,
        value="Mean",
        c="Mean",
        s="frac",
        cmap=cmap,
        verbose=1,
        **plot_kws,
    )
    for ax in cm1.heatmap_axes.ravel():
        ax.grid(
            axis="both",
            which="minor",
            color="grey",
            linestyle="--",
            alpha=0.6,
            zorder=0,
        )
    if outname is None:
        outname = f"{title}.pdf"
    plt.savefig(os.path.expanduser(outname), bbox_inches="tight")
    plt.show()
    return plot_data, df_cols, cm1


def get_genes_mean_frac(
    adata,
    obs=None,
    groupby="Subclass",
    modality="RNA",
    layer="mean",
    use_raw=False,
    expression_cutoff="p5",
    genes=None,
    normalize_per_cell=True,
    clip_norm_value=10,
    hypo_score=False,
):
    """Compute per-group mean expression and expressing-cell fraction.

    For single-cell data the mean and fraction are computed on the fly;
    for pseudobulk data with a pre-computed ``'mean'`` layer the values
    are read directly.

    Parameters
    ----------
    adata : AnnData or str
            Annotated data matrix or path to an ``.h5ad`` file.
    obs : str, DataFrame, or None, optional
            Cell metadata. Path to a TSV/CSV file or a DataFrame. If ``None``,
            uses ``adata.obs``.
    groupby : str, optional
            Column in *obs* to group cells by. Default is ``'Subclass'``.
    modality : str, optional
            Data modality (``'RNA'``, ``'ATAC'``, or methylation).
            Default is ``'RNA'``.
    layer : str, optional
            Layer name used to read pre-computed mean values from pseudobulk
            adata. Default is ``'mean'``.
    use_raw : bool, optional
            Whether to use ``adata.raw``. Default is ``False``.
    expression_cutoff : str or float, optional
            Threshold for a gene to be considered "expressed". Can be a
            percentile string (``'p5'``), ``'median'``, ``'mean'``, or a
            numeric value. Default is ``'p5'``.
    genes : list of str
            Gene names to compute statistics for.
    normalize_per_cell : bool, optional
            Whether to normalise methylation per cell. Default is ``True``.
    clip_norm_value : float, optional
            Clipping value for per-cell normalisation. Default is ``10``.
    hypo_score : bool, optional
            Whether to compute hypomethylation score. Default is ``False``.

    Returns
    -------
    pandas.DataFrame
            Tidy DataFrame with columns ``[groupby, 'Gene', 'Mean', 'frac']``.
    """
    assert genes is not None, "Please provide genes to plot."
    # adata could be single cell level or pseudobulk level (adata.layers['frac'] should be existed)
    if isinstance(adata, str):
        adata = anndata.read_h5ad(os.path.expanduser(adata), backed="r")
    all_vars = set(adata.var_names.tolist())
    keep_genes = list(
        set(all_vars) & set(genes)
    )  # keep_genes=[g for g in all_vars if g in genes]
    error_genes = [g for g in genes if g not in keep_genes]
    if len(error_genes) > 0:
        logger.debug(f"genes not found in adata: {error_genes}")
    use_adata = adata[:, keep_genes].to_memory()  # type: ignore
    if adata.isbacked:
        adata.file.close()  # close the file to save memory
    if "mean" not in use_adata.layers:  # raw count of single cell level adata
        # calculate mean and frac for each gene from single cell data
        if use_raw and use_adata.raw is not None:
            # use_adata.X=use_adata.raw.X.copy()
            use_adata_raw = use_adata.raw[:, use_adata.var_names.tolist()].to_adata()
            use_adata.X = use_adata_raw[
                use_adata.obs_names.tolist(), use_adata.var_names.tolist()
            ].X.copy()  # type: ignore
            del use_adata_raw
        if obs is not None:
            if isinstance(obs, str):
                sep = "\t" if obs.endswith(".tsv") or obs.endswith(".txt") else ","
                obs = pd.read_csv(os.path.expanduser(obs), sep=sep, index_col=0)
            assert isinstance(obs, pd.DataFrame), (
                "obs should be a dataframe or a path to a csv/tsv file."
            )
        else:
            obs = use_adata.obs.copy()
        overlapped_cells = list(
            set(use_adata.obs_names.tolist()) & set(obs.index.tolist())
        )
        obs = obs.loc[overlapped_cells]
        use_adata = use_adata[overlapped_cells, :]  # type: ignore

        if modality not in ["RNA", "ATAC"] and normalize_per_cell:
            use_adata = normalize_mc_by_cell(
                use_adata=use_adata,
                normalize_per_cell=normalize_per_cell,
                clip_norm_value=clip_norm_value,
                hypo_score=hypo_score,
            )

        data = use_adata.to_df()  # rows are cells or cell types, columns are genes
        if modality in ["RNA", "ATAC"] and isinstance(expression_cutoff, str):
            if expression_cutoff == "median":
                cutoff = data.stack().median()
            elif expression_cutoff == "mean":
                cutoff = data.stack().mean()
            else:  # quantile, such as p5,p95
                f = float(expression_cutoff.replace("p", ""))
                cutoff = data.stack().quantile(f / 100)
            expression_cutoff = cutoff

        data[groupby] = obs.loc[data.index.tolist(), groupby].tolist()  # type: ignore
        plot_data = data.groupby(groupby).mean().stack().reset_index()
        plot_data.columns = [groupby, "Gene", "Mean"]
        if "frac" in use_adata.layers:
            D = use_adata.to_df(layer="frac").stack().to_dict()
        else:
            if modality not in ["RNA", "ATAC"]:  # methylation, cutoff = 1
                assert normalize_per_cell, "Normalized methylation fraction is required"
                hypo_frac = data.groupby(
                    groupby
                ).agg(
                    lambda x: x[x < 1].shape[0] / x.shape[0]
                )  # fraction of cells showing hypomethylation for the corresponding genes
                D = hypo_frac.stack().to_dict()
            else:  # for RNA
                logger.info(f"Using expression cutoff: {expression_cutoff}")
                frac = data.groupby(groupby).agg(
                    lambda x: x[x > expression_cutoff].shape[0] / x.shape[0]
                )  # raw count > expression_cutoff means the gene is expressed
                D = frac.stack().to_dict()
        plot_data["frac"] = (
            plot_data.loc[:, [groupby, "Gene"]]
            .apply(lambda x: tuple(x.tolist()), axis=1)
            .map(D)
        )
    else:
        plot_data = use_adata.to_df(layer=layer).stack().reset_index()
        plot_data.columns = [groupby, "Gene", "Mean"]
        D = use_adata.to_df(layer="frac").stack().to_dict()
        plot_data["frac"] = (
            plot_data.loc[:, [groupby, "Gene"]]
            .apply(lambda x: tuple(x.tolist()), axis=1)
            .map(D)
        )
    return plot_data


def interactive_dotHeatmap(
    adata=None,
    obs=None,
    genes=None,
    groupby="Subclass",
    modality="RNA",
    title=None,
    use_raw=False,
    expression_cutoff="p5",
    normalize_per_cell=True,
    clip_norm_value=10,
    width=900,
    height=700,
    gene_order=None,
    colorscale="greens",
    vmin="p1",
    vmax="p99",
    show=True,
    reversescale=False,
    size_min=5,
    size_max=30,
    renderer="notebook",
):
    """Create an interactive Plotly dot-heatmap for multiple genes.

    Dot size encodes the fraction of expressing cells and colour encodes
    the mean expression value.

    Parameters
    ----------
    adata : AnnData, str, or None, optional
            Annotated data matrix or path to an ``.h5ad`` file.
    obs : str, DataFrame, or None, optional
            Cell metadata.
    genes : list of str
            Genes to plot.
    groupby : str, optional
            Column in *obs* to group cells by. Default is ``'Subclass'``.
    modality : str, optional
            Data modality. Default is ``'RNA'``.
    title : str or None, optional
            Plot title. Defaults to *groupby*.
    use_raw : bool, optional
            Whether to use ``adata.raw``. Default is ``False``.
    expression_cutoff : str or float, optional
            Threshold for a gene to be considered expressed.
            Default is ``'p5'``.
    normalize_per_cell : bool, optional
            Whether to normalise methylation per cell. Default is ``True``.
    clip_norm_value : float, optional
            Clipping value for per-cell normalisation. Default is ``10``.
    width : int, optional
            Figure width in pixels. Default is ``900``.
    height : int, optional
            Figure height in pixels. Default is ``700``.
    gene_order : list of str or None, optional
            Explicit ordering of genes on the Y axis.
    colorscale : str, optional
            Plotly colour scale name. Default is ``'greens'``.
    vmin : str, optional
            Colour-range lower bound as a percentile string (e.g. ``'p1'``).
            Default is ``'p1'``.
    vmax : str, optional
            Colour-range upper bound as a percentile string (e.g. ``'p99'``).
            Default is ``'p99'``.
    show : bool, optional
            Whether to display the figure immediately. Default is ``True``.
    reversescale : bool, optional
            Whether to reverse the colour scale. Default is ``False``.
    size_min : int, optional
            Minimum marker size in pixels. Default is ``5``.
    size_max : int, optional
            Maximum marker size in pixels. Default is ``30``.
    renderer : str, optional
            Plotly renderer backend. Default is ``'notebook'``.

    Returns
    -------
    plotly.graph_objects.Figure or None
            The figure object if ``show=False``, otherwise ``None``.
    """
    import plotly.io as pio
    import plotly.graph_objects as go

    if renderer is not None:
        pio.renderers.default = renderer
    plot_data = get_genes_mean_frac(
        adata,
        obs=obs,
        groupby=groupby,
        modality=modality,
        use_raw=use_raw,
        expression_cutoff=expression_cutoff,
        genes=genes,
        normalize_per_cell=normalize_per_cell,
        clip_norm_value=clip_norm_value,
        hypo_score=False,
    )  # columns: [groupby,'Gene','Mean','frac']
    # Build a Plotly dot-heatmap using scatter markers on categorical axes.
    # x: groups (columns), y: genes (rows)
    x_labels = plot_data[groupby].unique().tolist()
    if gene_order is None:
        y_labels = plot_data["Gene"].unique().tolist()
    else:
        y_labels = [g for g in gene_order if g in plot_data["Gene"].unique()]

    # Ensure ordering
    plot_data["x_cat"] = pd.Categorical(plot_data[groupby], categories=x_labels)
    plot_data["y_cat"] = pd.Categorical(plot_data["Gene"], categories=y_labels)

    # marker sizes: scale 'frac' (0-1) to reasonable pixel sizes
    frac_vals = plot_data["frac"].fillna(0).astype(float)
    sizes = (frac_vals * (size_max - size_min) + size_min).tolist()

    # marker colors: use Mean
    mean_vals = plot_data["Mean"].astype(float).tolist()

    hover_text = [
        f"Group: {g}<br>Gene: {ge}<br>Mean: {m:.4g}<br>Frac: {f:.3g}"
        for g, ge, m, f in zip(
            plot_data[groupby].tolist(),
            plot_data["Gene"].tolist(),
            mean_vals,
            frac_vals,
        )
    ]
    vmin_quantile = float(int(vmin.replace("p", "")) / 100)
    vmax_quantile = float(int(vmax.replace("p", "")) / 100)
    marker_dict = dict(
        size=sizes,
        color=mean_vals,
        colorscale=colorscale,
        showscale=True,
        colorbar=dict(title="Mean"),
        reversescale=reversescale,
        sizemode="area",
        opacity=0.9,
        cmin=plot_data["Mean"].quantile(vmin_quantile),
        cmax=plot_data["Mean"].quantile(vmax_quantile),
    )

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=plot_data[groupby].tolist(),
            y=plot_data["Gene"].tolist(),
            mode="markers",
            marker=marker_dict,
            text=hover_text,
            hoverinfo="text",
        )
    )

    # Layout: categorical axes with explicit ordering
    fig.update_xaxes(
        type="category", categoryorder="array", categoryarray=x_labels, tickangle=-45
    )
    fig.update_yaxes(
        type="category", categoryorder="array", categoryarray=list(reversed(y_labels))
    )
    if title is None:
        title = groupby
    fig.update_layout(
        title=title or "",
        xaxis_title=groupby,
        yaxis_title="Gene",
        width=width,
        height=height,
        plot_bgcolor="white",
    )

    if show:
        filename = f"dotHeatmap.{groupby}"
        show_fig(fig, filename=filename)
    else:
        return fig


def get_boxplot_data(adata, variable, gene, obs=None):
    """Extract a two-column DataFrame of group labels and gene values.

    Parameters
    ----------
    adata : AnnData
            Annotated data matrix (may be backed).
    variable : str
            Column in *obs* that defines the grouping variable.
    gene : str
            Gene name whose expression values are extracted.
    obs : str, DataFrame, or None, optional
            Cell metadata. Path to a TSV/CSV file or a DataFrame.
            If ``None``, uses ``adata.obs``.

    Returns
    -------
    pandas.DataFrame
            DataFrame with columns ``[variable, gene]`` for overlapping cells.
    """
    assert isinstance(adata, anndata.AnnData)
    if adata.isbacked:  # type: ignore
        use_adata = adata[:, gene].to_memory()  # type: ignore
    else:
        use_adata = adata[:, gene].copy()  # type: ignore
    if isinstance(obs, str):
        obs_path = os.path.abspath(os.path.expanduser(obs))
        sep = "\t" if obs_path.endswith(".tsv") or obs_path.endswith(".txt") else ","
        data = pd.read_csv(obs_path, index_col=0, sep=sep)
    else:
        assert isinstance(obs, pd.DataFrame)
        data = obs.copy()
    overlap_idx = data.index.intersection(use_adata.obs_names)
    data = data.loc[overlap_idx]
    use_adata = use_adata[overlap_idx, :]  # type: ignore

    if gene is not None:
        data[gene] = use_adata.to_df()[gene].tolist()  # type: ignore
    return data.loc[:, [variable, gene]]


def has_stats(adata):
    """Check whether an AnnData object contains pre-computed summary statistics.

    Looks for layers ``'min'``, ``'q25'``, ``'q50'``, ``'q75'``, ``'max'``,
    ``'mean'``, and ``'std'``.

    Parameters
    ----------
    adata : AnnData or str
            Annotated data matrix or path to an ``.h5ad`` file.

    Returns
    -------
    bool
            ``True`` if all required stat layers are present.
    """
    if isinstance(adata, str):
        adata = anndata.read_h5ad(adata, backed="r")
    flag = True
    for k in ["min", "q25", "q50", "q75", "max", "mean", "std"]:
        if k not in adata.layers:
            flag = False
            break
    return flag


def plot_interactive_boxlot_from_data(
    adata,
    obs,
    variable,
    gene,
    palette=None,
    width=1100,
    height=700,
    title=None,
):
    """Create a Plotly boxplot from raw single-cell expression data.

    Parameters
    ----------
    adata : AnnData
            Annotated data matrix.
    obs : str or DataFrame
            Cell metadata used to extract group labels and gene expression.
    variable : str
            Grouping column in *obs*.
    gene : str
            Gene name to plot on the Y axis.
    palette : str, dict, or None, optional
            Color palette for box colours.
    width : int, optional
            Figure width in pixels. Default is ``1100``.
    height : int, optional
            Figure height in pixels. Default is ``700``.
    title : str or None, optional
            Plot title. Defaults to ``'Boxplot: {gene} by {variable}'``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.express as px

    plot_df = get_boxplot_data(adata, variable, gene, obs=obs)
    # Preserve existing Y-axis extreme filtering logic (remove 1% and 99% extremes)
    range_y = [plot_df[gene].quantile(0.01), plot_df[gene].quantile(0.99)]
    # color_discrete_map=get_colors(adata,variable,palette=palette)
    color_discrete_map = load_color_palette(
        palette=palette, adata=adata, groups=variable
    )
    if color_discrete_map is not None:
        keys = list(color_discrete_map.keys())  # type: ignore
        for k in keys:
            if k not in plot_df[variable].unique().tolist():
                del color_discrete_map[k]  # type: ignore
    if title is None:
        title = f"Boxplot: {gene} by {variable}"
    fig = px.box(
        plot_df,
        x=variable,
        y=gene,
        color=variable,
        color_discrete_sequence=px.colors.qualitative.D3,  # color palette (professional, unobtrusive)
        color_discrete_map=color_discrete_map,
        range_y=range_y,
        points=False,
        title=title,
        template="plotly_white",  # keep white background style
    )
    fig.update_xaxes(tickangle=-90, automargin=True)

    fig.update_traces(
        line_width=1.2,  # thinner lines for a more refined look
        notched=False,  # no notch, standard boxplot style
    )

    fig.update_layout(
        xaxis_title=variable,
        yaxis_title=gene,
        legend_title=variable,
        width=width,
        height=height,
    )
    return fig


def plot_interacrive_boxplot_from_stats(
    adata, variable, gene, palette=None, title=None, width=1100, height=700
):
    """Create a Plotly boxplot from pre-computed pseudobulk summary statistics.

    Constructs boxes using the ``'min'``, ``'q25'``, ``'q50'``, ``'q75'``,
    and ``'max'`` layers stored in *adata*.

    Parameters
    ----------
    adata : AnnData
            Pseudobulk AnnData with pre-computed stat layers.
    variable : str
            Grouping variable (used as X axis label).
    gene : str
            Gene name to plot.
    palette : str, dict, or None, optional
            Color palette for box colours.
    title : str or None, optional
            Plot title. Defaults to ``'Boxplot: {gene} by {variable}'``.
    width : int, optional
            Figure width in pixels. Default is ``1100``.
    height : int, optional
            Figure height in pixels. Default is ``700``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.express as px
    import plotly.graph_objects as go

    assert isinstance(adata, anndata.AnnData)
    if adata.isbacked:  # type: ignore
        use_adata = adata[:, gene].to_memory()  # type: ignore
    else:
        use_adata = adata[:, gene].copy()  # type: ignore
    if adata.isbacked:  # type: ignore
        adata.file.close()  # type: ignore

    stat_keys = ["min", "q25", "q50", "q75", "max", "mean", "std"]
    plot_data = []
    for k in stat_keys:
        df = use_adata.to_df(layer=k)[gene]
        df.name = k
        plot_data.append(df)
    plot_data = pd.concat(
        plot_data, axis=1
    )  # rows are cell types and columns are stats
    groups = plot_data.sort_values("q50").index.tolist()
    plot_data = plot_data.loc[groups]

    # build figure with one Box per group using precomputed quartiles/fences
    fig = go.Figure()
    # optional color mapping
    # color_discrete_map=get_colors(adata,variable,palette=palette)
    color_discrete_map = load_color_palette(
        palette=palette, adata=adata, groups=variable
    )
    palette = px.colors.qualitative.D3
    color = None
    for group, row in plot_data.iterrows():
        i = groups.index(group)
        q1 = row["q25"]
        med = row["q50"]
        q3 = row["q75"]
        low = row["min"]
        high = row["max"]
        if color_discrete_map is not None and group in color_discrete_map:
            color = color_discrete_map[group]
        else:
            color = palette[i % len(palette)]
        # Box from precomputed stats (single-element arrays)
        fig.add_trace(
            go.Box(
                x=[group],
                q1=[q1],
                median=[med],
                q3=[q3],
                lowerfence=[low],
                upperfence=[high],
                boxpoints=False,
                marker=dict(color=color),
                name=str(group),
                showlegend=True,
            )
        )
        # mean as a scatter point with std error bar
        # fig.add_trace(
        # 	go.Scatter(
        # 		x=[group],
        # 		y=[mean],
        # 		mode='markers',
        # 		marker=dict(symbol='diamond', size=8, color='black'),
        # 		error_y=dict(type='data', array=[std], visible=True),
        # 		name='mean',
        # 		showlegend=False
        # 	)
        # )
    if title is None:
        title = f"Boxplot: {gene} by {variable}"
    fig.update_xaxes(tickangle=-90, automargin=True)
    fig.update_layout(
        title=title,
        xaxis_title=variable,
        yaxis_title=gene,
        legend_title=variable,
        template="plotly_white",
        width=width,
        height=height,
    )
    return fig


def interactive_boxplot(
    adata,
    variable,
    gene,
    obs=None,
    palette=None,
    title=None,
    width=1100,
    height=700,
    show=True,
    renderer="notebook",
):
    """Create an interactive Plotly boxplot for a gene grouped by a variable.

    Automatically detects whether *adata* contains pre-computed summary
    statistics (pseudobulk) or raw data and dispatches to the appropriate
    plotting function.

    Parameters
    ----------
    adata : AnnData or str
            Annotated data matrix or path to an ``.h5ad`` file.
    variable : str
            Grouping column.
    gene : str
            Gene name to plot on the Y axis.
    obs : str, DataFrame, or None, optional
            Cell metadata. If ``None``, uses ``adata.obs``.
    palette : str, dict, or None, optional
            Color palette for box colours.
    title : str or None, optional
            Plot title.
    width : int, optional
            Figure width in pixels. Default is ``1100``.
    height : int, optional
            Figure height in pixels. Default is ``700``.
    show : bool, optional
            Whether to display the figure. Default is ``True``.
    renderer : str, optional
            Plotly renderer backend. Default is ``'notebook'``.

    Returns
    -------
    plotly.graph_objects.Figure or None
            The figure object if ``show=False``, otherwise ``None``.
    """
    import plotly.io as pio

    if renderer is not None:
        pio.renderers.default = renderer
    if isinstance(adata, str):
        adata = anndata.read_h5ad(adata, backed="r")
    else:
        assert isinstance(adata, anndata.AnnData)
    if obs is None:
        obs = adata.obs.copy()  # type: ignore
    if not has_stats(adata):
        fig = plot_interactive_boxlot_from_data(
            adata,
            obs,
            variable,
            gene,
            palette=palette,
            title=title,
            width=width,
            height=height,
        )
    else:  # pseudobulk level with precomputed stats
        fig = plot_interacrive_boxplot_from_stats(
            adata,
            variable,
            gene,
            palette=palette,
            title=title,
            width=width,
            height=height,
        )
    if show:
        filename = f"boxplot.{variable}.{gene}"
        show_fig(fig, filename=filename)
        return None
    else:
        return fig
