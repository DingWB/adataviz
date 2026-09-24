"""Embedding scatter plots: ``plot_categorical`` and ``plot_gene``.

These wrap the lower-level helpers in :mod:`adataviz.utils` (the
``categorical_scatter`` / ``continuous_scatter`` engines) so that they
can be used directly from ``adataviz.pl`` with palette resolution from
``adata.uns`` and a consistent ``legend_kws`` interface.

``plot_gene`` is now scoped to the embedding plot only; for box / strip
/ violin plots of single genes use :func:`adataviz.pl.boxplot` or
:func:`adataviz.pl.stacked_violinplot`.
"""

from __future__ import annotations

import os
from typing import Any, Mapping, Optional

import anndata
import matplotlib.pyplot as plt
import pandas as pd

from ..tools import load_adata
from ..utils import (
    categorical_scatter,
    continuous_scatter,
    normalize_mc_by_cell,
)

__all__ = ["plot_categorical", "plot_gene", "assign_default_colors"]


# ---------------------------------------------------------------------------
# Default categorical palette
# ---------------------------------------------------------------------------


_DEFAULT_20 = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
    "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
    "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
]
_DEFAULT_28 = _DEFAULT_20 + [
    "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#bd9e39",
    "#ad494a",
]


def assign_default_colors(adata, key: str) -> dict:
    """Assign default categorical colours for ``adata.obs[key]``.

    Mirrors scanpy's palette logic: the current matplotlib colour cycle
    (``rcParams["axes.prop_cycle"]``) is used when it has enough colours,
    otherwise a built-in 20-colour or 28-colour palette is used, falling
    back to a cycled ``tab20`` colormap when there are more than 28
    categories. The resulting list is written to
    ``adata.uns[f"{key}_colors"]`` (aligned with ``cat.categories``) so it
    can be reused by downstream plotting functions.

    Parameters
    ----------
    adata : AnnData
        Annotated data object. ``adata.obs[key]`` must be categorical (or
        expose a ``.cat.categories`` accessor).
    key : str
        Name of the categorical column in ``adata.obs`` to assign colours
        for.

    Returns
    -------
    dict
        Mapping of ``{category: hex_colour}`` for every category in
        ``adata.obs[key]``. The same colours (as an ordered list) are also
        stored in ``adata.uns[f"{key}_colors"]`` as a side effect.
    """
    from matplotlib.colors import to_hex
    from matplotlib import rcParams

    categories = adata.obs[key].cat.categories
    n = len(categories)
    cycle = rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(cycle) >= n:
        palette = [to_hex(c) for c in cycle[:n]]
    elif n <= 20:
        palette = _DEFAULT_20[:n]
    elif n <= 28:
        palette = _DEFAULT_28[:n]
    else:
        # Fall back to tab20 cycled
        cmap = plt.get_cmap("tab20")
        palette = [to_hex(cmap(i % cmap.N)) for i in range(n)]
    adata.uns[f"{key}_colors"] = palette
    return {cat: col for cat, col in zip(categories, palette)}


# ---------------------------------------------------------------------------
# plot_categorical
# ---------------------------------------------------------------------------


def plot_categorical(
    adata,
    ax=None,
    basis: str = "umap",
    groupby: str = "Subclass",
    coding: bool = True,
    coded_marker: bool = True,
    save: Optional[str] = None,
    palette=None,
    sheet_name: Optional[str] = None,
    show: bool = True,
    figsize=(4, 3.5),
    dpi: int = 300,
    pad_inches: float = 0.02,
    ncol: Optional[int] = None,
    fontsize: int = 5,
    legend_fontsize: int = 5,
    legend_kws: Optional[Mapping[str, Any]] = None,
    legend_title_fontsize: int = 5,
    marker_fontsize: int = 4,
    marker_pad: float = 0.1,
    linewidth: float = 0.5,
    axis_format: str = "tiny",
    alpha: float = 0.7,
    text_kws: Optional[Mapping[str, Any]] = None,
    **kwargs,
):
    """Categorical scatter on an embedding (e.g. UMAP) coloured by ``groupby``.

    Colours are resolved in the following order of precedence: the
    ``palette`` argument (dict or Excel path), an existing
    ``adata.uns[f"{groupby}_colors"]`` entry, or a freshly generated
    default palette via :func:`assign_default_colors`. The resolved
    colours are cached back into ``adata.uns`` and the plot is drawn with
    the :func:`adataviz.utils.categorical_scatter` engine.

    Parameters
    ----------
    adata : AnnData or str
        Annotated data object, or a path passed to
        :func:`adataviz.tools.load_adata`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If ``None`` (default) a new figure and axes are
        created using ``figsize`` and ``dpi``.
    basis : str, default ``"umap"``
        Embedding key in ``adata.obsm``. A leading ``X_`` prefix is
        stripped automatically (e.g. ``"X_umap"`` becomes ``"umap"``).
    groupby : str, default ``"Subclass"``
        Name of the categorical column in ``adata.obs`` used to colour and
        annotate the points. Non-categorical columns are coerced to
        ``category`` dtype.
    coding : bool, default ``True``
        Forwarded to the scatter engine as ``coding``; when ``True`` a
        short numeric code is used in place of the full category label for
        the on-plot text annotation.
    coded_marker : bool, default ``True``
        Forwarded to the scatter engine as ``coded_marker``; controls
        whether the numeric code marker is drawn next to each cluster when
        coded labels are in use.
    save : str, optional
        Path to save the figure to (``~`` is expanded). Saved with
        ``bbox_inches="tight"`` and the given ``pad_inches``. If ``None``
        (default) the figure is not saved.
    palette : dict or str, optional
        Colour source. A ``dict`` maps ``{category: colour}`` directly. A
        ``str`` is treated as a path to an Excel file read with
        ``sheet_name`` and an index column, using its ``Hex`` column;
        categories absent from the data are dropped and missing ones
        default to ``"gray"``. If ``None`` (default) colours fall back to
        ``adata.uns`` or a generated palette.
    sheet_name : str, optional
        Worksheet name used when ``palette`` is an Excel path. Defaults to
        ``groupby`` when ``None``.
    show : bool, default ``True``
        Whether to call ``plt.show()`` after drawing.
    figsize : tuple of float, default ``(4, 3.5)``
        Figure size in inches, used only when ``ax`` is ``None``.
    dpi : int, default ``300``
        Resolution of the created figure and of the saved output.
    pad_inches : float, default ``0.02``
        Padding (in inches) around the figure when saving with
        ``bbox_inches="tight"``.
    ncol : int, optional
        Number of legend columns. If ``None`` (default) it is computed
        automatically as roughly one column per 20 categories so large
        legends do not overflow the figure.
    fontsize : int, default ``5``
        Default font size for on-plot text annotations (used as the
        default for ``text_kws["fontsize"]``).
    legend_fontsize : int, default ``5``
        Font size of the legend entry labels.
    legend_kws : Mapping, optional
        Extra keyword arguments forwarded to the legend. Sensible
        defaults are filled in for ``fontsize``, ``title`` (set to
        ``groupby``), ``title_fontsize`` and ``ncol``.
    legend_title_fontsize : int, default ``5``
        Font size of the legend title.
    marker_fontsize : int, default ``4``
        Font size of the coded cluster markers (forwarded to the engine).
    marker_pad : float, default ``0.1``
        Padding around coded cluster markers (forwarded to the engine).
    linewidth : float, default ``0.5``
        Line width of the scatter marker edges (forwarded to the engine).
    axis_format : str, default ``"tiny"``
        Axis styling passed to the engine; ``"tiny"`` draws small
        embedding-style axis arrows/labels instead of full axes.
    alpha : float, default ``0.7``
        Opacity of the scatter points (forwarded to the engine).
    text_kws : Mapping, optional
        Extra keyword arguments forwarded to the text annotation drawing.
        The ``fontsize`` key defaults to ``fontsize`` when not provided.
    **kwargs
        Additional keyword arguments forwarded to
        :func:`adataviz.utils.categorical_scatter`. Several defaults are
        pre-populated (``hue``, ``text_anno``, ``text_kws``, ``luminance``,
        ``dodge_text``, ``axis_format``, ``show_legend``,
        ``marker_fontsize``, ``marker_pad``, ``linewidth``, ``alpha`` and
        ``dodge_kws``); ``coding`` and ``coded_marker`` are set
        unconditionally. Useful pass-through options include ``s``
        (marker size), ``sizes``/``size``/``size_norm`` (size mapping),
        ``zoomxy`` (axis zoom factor), ``max_points`` (subsampling cap),
        ``rasterized`` and ``outline``.

    Returns
    -------
    matplotlib.axes.Axes
        The axes the categorical scatter was drawn on.
    """
    if basis.startswith("X_"):
        basis = basis[2:]
    if sheet_name is None:
        sheet_name = groupby
    adata = load_adata(adata)
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
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
                existed = adata.obs[groupby].unique().tolist()
                for k in list(colors.keys()):
                    if k not in existed:
                        del colors[k]
                for k in existed:
                    colors.setdefault(k, "gray")
            except Exception:
                colors = None
    if colors is None:
        if f"{groupby}_colors" in adata.uns:
            colors = {
                cat: col
                for cat, col in zip(
                    adata.obs[groupby].cat.categories.tolist(),
                    adata.uns[f"{groupby}_colors"],
                )
            }
        else:
            colors = assign_default_colors(adata, groupby)
    else:
        adata.uns[f"{groupby}_colors"] = [
            colors.get(k, "grey") for k in adata.obs[groupby].cat.categories.tolist()
        ]

    text_kws = {} if text_kws is None else dict(text_kws)
    text_kws.setdefault("fontsize", fontsize)
    kwargs.setdefault("hue", groupby)
    kwargs.setdefault("text_anno", groupby)
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

    legend_kws = {} if legend_kws is None else dict(legend_kws)
    legend_kws.setdefault("fontsize", legend_fontsize)
    legend_kws.setdefault("title", groupby)
    legend_kws.setdefault("title_fontsize", legend_title_fontsize)
    if ncol is None:
        # Wrap large categorical legends into multiple columns (~20 per
        # column) so they don't run off the bottom of the figure.
        n_cats = int(adata.obs[groupby].nunique(dropna=True))
        ncol = max(1, (n_cats + 19) // 20)
    legend_kws.setdefault("ncol", ncol)
    kwargs.setdefault(
        "dodge_kws",
        dict(
            arrowprops=dict(
                arrowstyle="->",
                fc="grey",
                ec="none",
                connectionstyle="angle,angleA=-90,angleB=180,rad=5",
            ),
            autoalign="xy",
        ),
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure
    categorical_scatter(
        data=adata[adata.obs[groupby].notna()],
        ax=ax,
        basis=basis,
        palette=colors,
        legend_kws=legend_kws,
        **kwargs,
    )
    if save:
        fig.set_dpi(dpi)
        fig.savefig(
            os.path.expanduser(save),
            dpi=dpi,
            bbox_inches="tight",
            pad_inches=pad_inches,
        )
    if show:
        plt.show()
    return ax


# ---------------------------------------------------------------------------
# plot_gene (UMAP only)
# ---------------------------------------------------------------------------


def plot_gene(
    adata,
    gene: str,
    basis: str = "umap",
    obs=None,
    query_str: Optional[str] = None,
    title: Optional[str] = None,
    cmap: str = "parula",
    hue_norm=None,
    cbar_kws=None,
    axis_format: str = "tiny",
    obsm=None,
    normalize_per_cell: bool = False,
    hypo_score: bool = False,
    clip_norm_value: float = 10,
    figsize=(4, 4),
    dpi: int = 300,
    pad_inches: float = 0.02,
    ax=None,
    save: Optional[str] = None,
    show: bool = True,
    scatter_kws: Optional[Mapping[str, Any]] = None,
):
    """Plot a single gene's expression on an embedding.

    The previous ``plot_gene`` produced both a scatter and a per-group
    violin/box plot. The violin/box behaviour now lives in
    :func:`adataviz.pl.boxplot`; this function is purely about the
    embedding visualisation. The single gene is sliced out (loading from
    disk in ``backed`` mode when ``adata`` is a path), optionally
    normalised and re-embedded, then drawn with the
    :func:`adataviz.utils.continuous_scatter` engine.

    Parameters
    ----------
    adata : AnnData or str
        Source data, or a path to an ``.h5ad`` file (opened in ``backed``
        mode and closed after the gene is extracted).
    gene : str
        Gene name; must exist in ``adata.var_names``. Its expression is
        added to ``obs`` and used as the colour ``hue``.
    basis : str, default ``"umap"``
        Embedding key in ``adata.obsm``. A leading ``X_`` prefix is
        stripped automatically.
    obs : pandas.DataFrame or str, optional
        Replacement observation metadata, either a DataFrame or a path to
        a tab-separated file (index in the first column). When ``None``
        (default) the object's own ``obs`` is used. Only cells whose
        indices overlap the expression data are kept.
    query_str : str, optional
        Pandas ``query`` expression applied to ``obs`` to subset cells
        before plotting. If ``None`` (default) no filtering is applied.
    title : str, optional
        Axes title. If ``None`` (default) no title is set.
    cmap : str, default ``"parula"``
        Colormap name forwarded to the continuous scatter engine.
    hue_norm : optional
        Normalisation for the colour scale (e.g. a
        ``matplotlib.colors.Normalize`` instance or ``(vmin, vmax)``
        tuple) forwarded to the engine.
    cbar_kws : dict, optional
        Keyword arguments for the colorbar. Defaults to
        ``dict(extendfrac=0.1)`` when ``None``.
    axis_format : str, default ``"tiny"``
        Axis styling passed to the engine; ``"tiny"`` draws small
        embedding-style axes.
    obsm : AnnData or str, optional
        Alternative embedding source, either an AnnData or a path to an
        ``.h5ad`` file. When provided, its ``obsm`` (and any extra ``obs``
        columns) replace those of the sliced data for cells present in
        both.
    normalize_per_cell : bool, default ``False``
        Whether to normalise methylation/expression per cell via
        :func:`adataviz.utils.normalize_mc_by_cell` before plotting.
    hypo_score : bool, default ``False``
        Passed to the normalisation step; when ``True`` values are turned
        into a hypo-methylation score.
    clip_norm_value : float, default ``10``
        Upper clipping value applied during per-cell normalisation.
    figsize : tuple of float, default ``(4, 4)``
        Figure size in inches, used only when ``ax`` is ``None``.
    dpi : int, default ``300``
        Resolution of the created figure and of the saved output.
    pad_inches : float, default ``0.02``
        Padding (in inches) around the figure when saving with
        ``bbox_inches="tight"``.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If ``None`` (default) a new figure and axes are
        created using ``figsize`` and ``dpi``.
    save : str, optional
        Path to save the figure to (``~`` is expanded). Saved with
        ``bbox_inches="tight"`` and the given ``pad_inches``. If ``None``
        (default) the figure is not saved.
    show : bool, default ``True``
        Whether to call ``plt.show()`` after drawing.
    scatter_kws : Mapping, optional
        Extra keyword arguments forwarded to
        :func:`adataviz.utils.continuous_scatter`. Defaults to an empty
        dict. Useful pass-through options include ``s`` (marker size),
        ``alpha`` (via ``scatter_kws``), ``zoomxy`` (axis zoom factor),
        ``max_points`` (subsampling cap), ``rasterized``, ``size``/
        ``size_norm``/``sizes`` (size mapping) and ``outline``.

    Returns
    -------
    matplotlib.axes.Axes
        The axes the gene expression scatter was drawn on.
    """
    if cbar_kws is None:
        cbar_kws = dict(extendfrac=0.1)
    if scatter_kws is None:
        scatter_kws = {}

    raw = anndata.read_h5ad(os.path.expanduser(adata), backed="r") if isinstance(
        adata, str
    ) else adata
    use = raw[:, gene].to_memory() if raw.isbacked else raw[:, gene].copy()
    if isinstance(adata, str) and raw.isbacked:
        raw.file.close()

    if normalize_per_cell:
        use = normalize_mc_by_cell(
            use_adata=use,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
            hypo_score=hypo_score,
        )

    if obsm is not None:
        if isinstance(obsm, str):
            obsm = anndata.read_h5ad(os.path.expanduser(obsm), backed="r")
        keep = list(set(use.obs_names) & set(obsm.obs_names))
        use = use[keep, :]
        use.obsm = obsm[keep].obsm
        for col in obsm.obs.columns:
            if col not in use.obs.columns:
                use.obs[col] = obsm.obs.loc[use.obs_names, col].tolist()
        if hasattr(obsm, "isbacked") and obsm.isbacked:
            obsm.file.close()

    if obs is not None:
        if isinstance(obs, str):
            obs = pd.read_csv(os.path.expanduser(obs), sep="\t", index_col=0)
        else:
            obs = obs.copy()
    else:
        obs = use.obs.copy()
    if query_str is not None:
        obs = obs.query(query_str)
    overlap = list(set(use.obs_names) & set(obs.index))
    use = use[overlap, :]
    use.obs = obs.loc[use.obs_names]
    use.obs[gene] = use.to_df().loc[use.obs_names, gene].tolist()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure
    continuous_scatter(
        data=use,
        ax=ax,
        cmap=cmap,
        hue_norm=hue_norm,
        cbar_kws=cbar_kws,
        hue=gene,
        axis_format=axis_format,
        text_anno=None,
        basis=basis if not basis.startswith("X_") else basis[2:],
        **scatter_kws,
    )
    if not title is None:
        ax.set_title(title)
    if save:
        fig.set_dpi(dpi)
        fig.savefig(
            os.path.expanduser(save),
            dpi=dpi,
            bbox_inches="tight",
            pad_inches=pad_inches,
        )
    if show:
        plt.show()
    return ax
