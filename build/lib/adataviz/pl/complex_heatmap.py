"""High-level wrappers around :mod:`PyComplexHeatmap`.

- :func:`complex_heatmap` -- aggregate ``adata`` cells by ``groupby``
  (mean expression of ``genes``), then draw a clustered/annotated heatmap.
- :func:`complex_dotplot` -- same aggregation, but rendered as a dot plot
  where dot size = fraction expressing and colour = mean expression.

Both wrappers expose:

- ``annotate`` -- list of extra ``adata.obs`` columns to add as row
  annotations (palettes auto-resolved from ``adata.uns``).
- ``plot_kws`` -- forwarded directly to the underlying
  ``ClusterMapPlotter`` / ``DotClustermapPlotter``.
- ``annot_kws`` -- forwarded to ``HeatmapAnnotation``.
- ``**kwargs`` -- final fall-through into the plotter, so any
  PyComplexHeatmap parameter not exposed explicitly is still reachable.

The wrappers also strip the spurious white outline that PyComplexHeatmap
sometimes leaves at the bottom of colorbars.
"""

from __future__ import annotations

import contextlib
import warnings
from typing import Any, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._utils import resolve_palette, save_or_show

__all__ = ["complex_heatmap", "complex_dotplot"]


@contextlib.contextmanager
def _quiet_pch(fig):
    """Silence upstream PyComplexHeatmap warnings while it draws ``fig``.

    - Disables ``figure.autolayout`` on the target figure so that
      PyComplexHeatmap's manual GridSpec layout does not trigger the
      "Axes are not compatible with tight_layout" warning on canvas draw.
    - Filters the pandas ``FutureWarning`` about callable ``mean`` /
      ``DataFrame.applymap`` that comes from ``PyComplexHeatmap.dotHeatmap``.
    """
    # Force this figure off any auto-layout engine (tight/constrained) that
    # ``figure.autolayout=True`` in our rcParams would otherwise install.
    try:
        fig.set_layout_engine("none")
    except Exception:
        pass
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=FutureWarning,
            module=r"PyComplexHeatmap\..*",
        )
        warnings.filterwarnings(
            "ignore",
            message=r".*not compatible with tight_layout.*",
            category=UserWarning,
        )
        yield


def _aggregate(adata, groupby, genes, layer=None, use_raw=False, expression_cutoff=0):
    import anndata as _ad

    if not isinstance(adata, _ad.AnnData):
        raise TypeError("complex_heatmap/complex_dotplot require AnnData input.")
    src = adata.raw if use_raw and adata.raw is not None else adata
    var_names = list(src.var_names)
    missing = [g for g in genes if g not in var_names]
    if missing:
        raise KeyError(f"genes not in adata.var_names: {missing[:5]}")
    idx = [var_names.index(g) for g in genes]
    if layer is not None and layer in adata.layers:
        X = adata.layers[layer][:, idx]
    else:
        X = src.X[:, idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)
    obs = adata.obs[[groupby]].copy()
    obs[groupby] = obs[groupby].astype(str)
    df = pd.DataFrame(X, index=obs.index, columns=list(genes))
    g = obs[groupby]
    mean_df = df.groupby(g, observed=True).mean()
    frac_df = (df > expression_cutoff).groupby(g, observed=True).mean()
    return mean_df, frac_df


def _strip_cbar_white_lines(cm):
    """Hide the residual white spines around PyComplexHeatmap colorbars."""
    for cbar in getattr(cm, "cbars", []) or []:
        try:
            cbar.outline.set_visible(False)
            for sp in cbar.ax.spines.values():
                sp.set_visible(False)
            cbar.ax.tick_params(length=2, width=0.5)
        except Exception:
            pass
    for cax in getattr(cm, "legend_axes", []) or []:
        try:
            for sp in cax.spines.values():
                sp.set_visible(False)
        except Exception:
            pass


def _row_annotation(
    adata, cats, groupby, palette, annotate, annotate_palettes, annot_kws
):
    import PyComplexHeatmap as pch

    primary_colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=adata, groupby=groupby
    )
    annos = {
        "Group": pch.anno_simple(
            pd.Series(cats, index=cats, name=groupby),
            colors=primary_colors,
            add_text=False,
            legend=False,
        )
    }
    if annotate:
        annotate_palettes = dict(annotate_palettes or {})
        sub = adata.obs[[groupby] + list(annotate)].copy()
        sub[groupby] = sub[groupby].astype(str)
        per_group = (
            sub.groupby(groupby, observed=True)
            .agg(lambda v: v.mode().iloc[0] if not v.mode().empty else v.iloc[0])
            .reindex(cats)
        )
        for col in annotate:
            if col not in per_group.columns:
                continue
            ann_cats = per_group[col].astype(str).tolist()
            uniq = list(dict.fromkeys(ann_cats))
            cpal = resolve_palette(
                annotate_palettes.get(col),
                uniq,
                sheet_name=col,
                adata=adata,
                groupby=col,
            )
            annos[col] = pch.anno_simple(
                pd.Series(ann_cats, index=cats, name=col),
                colors=cpal,
                add_text=False,
                legend=True,
            )
    annot_kws = dict(annot_kws or {})
    annot_kws.setdefault("axis", 0)
    annot_kws.setdefault("verbose", 0)
    return pch.HeatmapAnnotation(**annos, **annot_kws), primary_colors


def complex_heatmap(
    adata,
    genes: Sequence[str],
    groupby: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    z_score: Optional[str] = "row",
    cmap: str = "RdBu_r",
    palette: Union[None, Mapping[str, str], str] = None,
    annotate: Optional[Sequence[str]] = None,
    annotate_palettes: Optional[Mapping[str, Any]] = None,
    row_cluster: bool = True,
    col_cluster: bool = False,
    figsize=None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
    show_rownames: bool = True,
    show_colnames: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
    plot_kws: Optional[Mapping[str, Any]] = None,
    annot_kws: Optional[Mapping[str, Any]] = None,
    **kwargs,
):
    """Annotated clustered heatmap of mean gene expression per group.

    Cells in ``adata`` are aggregated by ``groupby`` (mean expression of the
    selected ``genes`` per group), optionally z-scored across rows or columns,
    and drawn as a clustered, annotated heatmap via
    :class:`PyComplexHeatmap.ClusterMapPlotter`. Each row is a group and each
    column a gene; a left-side ``Group`` annotation stripe (plus any extra
    ``annotate`` columns) labels the rows, and colour encodes the (standardised)
    mean expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must be a real :class:`~anndata.AnnData`; a
        ``TypeError`` is raised otherwise.
    genes : sequence of str
        Gene / feature names to display as heatmap columns. Every entry must
        exist in ``adata.var_names`` (or ``adata.raw.var_names`` when
        ``use_raw=True``); a ``KeyError`` is raised for missing genes.
    groupby : str
        Column in ``adata.obs`` used to aggregate cells into groups (heatmap
        rows). Values are coerced to ``str`` before grouping.
    layer : str, optional
        Name of an ``adata.layers`` entry to read expression from. If ``None``
        (default) or the name is absent, ``adata.X`` (or ``adata.raw.X`` when
        ``use_raw=True``) is used instead.
    use_raw : bool, default ``False``
        If ``True`` and ``adata.raw`` is present, read expression from
        ``adata.raw`` instead of ``adata``.
    z_score : {"row", "col", None}, default ``"row"``
        Standardisation applied to the aggregated mean matrix before plotting.
        ``"row"`` z-scores each group (subtract row mean, divide by row std),
        ``"col"`` z-scores each gene, and ``None`` leaves the raw means. Zero
        standard deviations are handled safely (result filled with 0).
    cmap : str, default ``"RdBu_r"``
        Matplotlib colormap name used for the heatmap fill.
    palette : dict, str or None, optional
        Colour spec for the primary ``Group`` (``groupby``) annotation stripe.
        Accepts an explicit ``{category: colour}`` mapping, a path to an Excel
        palette sheet, or ``None`` to auto-resolve from
        ``adata.uns[f"{groupby}_colors"]`` (falling back to matplotlib
        ``tab10``/``tab20``). See :func:`resolve_palette`.
    annotate : sequence of str, optional
        Extra ``adata.obs`` columns drawn as additional left-side row
        annotation stripes (one stripe per column). For each group the dominant
        (modal) value of the column is used.
    annotate_palettes : dict, optional
        Per-annotation-column palette spec, mapping an ``annotate`` column name
        to a palette (dict / Excel path / ``None``) resolved the same way as
        ``palette``.
    row_cluster : bool, default ``True``
        Whether to hierarchically cluster and reorder the rows (groups).
    col_cluster : bool, default ``False``
        Whether to hierarchically cluster and reorder the columns (genes).
    figsize : tuple of float, optional
        Figure size in inches. If ``None`` (default) it is derived from the
        number of genes and groups so tick labels stay readable on large
        panels.
    title : str, optional
        Figure suptitle drawn above the heatmap. No title if ``None``.
    save : str, optional
        Path to write the figure to (parent directories created as needed).
        If ``None`` (default) the figure is not saved.
    show : bool, default ``False``
        If ``True`` display the figure with ``plt.show()``.
    show_rownames : bool, default ``True``
        Whether to draw row (group) labels.
    show_colnames : bool, default ``True``
        Whether to draw column (gene) labels.
    legend_kws : dict, optional
        Extra keyword arguments forwarded as the ``legend_kws`` of the
        underlying ``ClusterMapPlotter`` (colorbar / legend styling).
    plot_kws : dict, optional
        Extra keyword arguments forwarded to
        :class:`PyComplexHeatmap.ClusterMapPlotter`, overriding the wrapper's
        defaults for any key they set.
    annot_kws : dict, optional
        Extra keyword arguments forwarded to
        :class:`PyComplexHeatmap.HeatmapAnnotation` used to build the row
        annotations.
    **kwargs
        Additional keyword arguments merged into ``plot_kws`` and passed
        straight through to ``ClusterMapPlotter``, so any PyComplexHeatmap
        parameter not exposed explicitly here remains reachable.

    Returns
    -------
    PyComplexHeatmap.ClusterMapPlotter
        The fitted plotter object holding the drawn heatmap (its ``fig`` /
        axes can be used for further customisation).
    """
    import PyComplexHeatmap as pch

    mean_df, _ = _aggregate(adata, groupby, list(genes), layer=layer, use_raw=use_raw)
    if z_score == "row":
        mean_df = (
            mean_df.sub(mean_df.mean(axis=1), axis=0)
            .div(mean_df.std(axis=1).replace(0, np.nan), axis=0)
            .fillna(0)
        )
    elif z_score == "col":
        mean_df = (
            mean_df.sub(mean_df.mean(axis=0), axis=1)
            .div(mean_df.std(axis=0).replace(0, np.nan), axis=1)
            .fillna(0)
        )
    cats = list(mean_df.index)
    row_anno, _ = _row_annotation(
        adata, cats, groupby, palette, annotate, annotate_palettes, annot_kws
    )

    # Scale the figure with the number of genes (columns) and groups (rows)
    # so tick labels stay readable instead of overlapping on large panels.
    if figsize is None:
        n_col, n_row = len(mean_df.columns), len(mean_df.index)
        figsize = (max(4.0, 0.28 * n_col + 2.5), max(3.0, 0.28 * n_row + 1.8))

    fig = plt.figure(figsize=figsize)
    plot_kws = dict(plot_kws or {})
    plot_kws.update(kwargs)
    # Right-side row labels (row_names_side="right") extend past the heatmap
    # into the legend column; widen the legend gap for long group names so
    # the colour legend does not overlap the y tick labels.
    _maxlen = max((len(str(c)) for c in cats), default=0)
    _auto_hpad = max(4, round(0.75 * _maxlen))
    base = dict(
        data=mean_df,
        cmap=cmap,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        show_rownames=show_rownames,
        show_colnames=show_colnames,
        left_annotation=row_anno,
        col_names_side="bottom",
        row_names_side="right",
        xticklabels_kws=dict(labelrotation=-45, labelsize=8, bottom=True),
        yticklabels_kws=dict(labelsize=8, right=True),
        legend_hpad=_auto_hpad,
        legend_side="right",
        legend_anchor="ax_heatmap",
        verbose=0,
    )
    for k, v in plot_kws.items():
        base[k] = v
    if legend_kws:
        base["legend_kws"] = dict(legend_kws)
    with _quiet_pch(fig):
        cm = pch.ClusterMapPlotter(**base)
    _strip_cbar_white_lines(cm)
    if title:
        fig.suptitle(title, y=1.02)
    with _quiet_pch(fig):
        save_or_show(fig, save, show=show)
    return cm


def complex_dotplot(
    adata,
    genes: Sequence[str],
    groupby: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    expression_cutoff: float = 0,
    cmap: str = "Reds",
    palette: Union[None, Mapping[str, str], str] = None,
    annotate: Optional[Sequence[str]] = None,
    annotate_palettes: Optional[Mapping[str, Any]] = None,
    row_cluster: bool = False,
    col_cluster: bool = False,
    figsize=None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
    show_rownames: bool = True,
    show_colnames: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
    plot_kws: Optional[Mapping[str, Any]] = None,
    annot_kws: Optional[Mapping[str, Any]] = None,
    marker: str = "o",
    **kwargs,
):
    """Dot heatmap: dot size = fraction expressing, colour = mean expression.

    Cells in ``adata`` are aggregated by ``groupby`` and the selected ``genes``
    are rendered as a dot plot via
    :class:`PyComplexHeatmap.DotClustermapPlotter`. Each row is a group and each
    column a gene; the dot **size** encodes the fraction of cells in the group
    expressing the gene (expression strictly above ``expression_cutoff``) and
    the dot **colour** encodes the group's mean expression. A left-side
    ``Group`` annotation stripe (plus any extra ``annotate`` columns) labels the
    rows.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must be a real :class:`~anndata.AnnData`; a
        ``TypeError`` is raised otherwise.
    genes : sequence of str
        Gene / feature names to display as dot-plot columns. Every entry must
        exist in ``adata.var_names`` (or ``adata.raw.var_names`` when
        ``use_raw=True``); a ``KeyError`` is raised for missing genes.
    groupby : str
        Column in ``adata.obs`` used to aggregate cells into groups (dot-plot
        rows). Values are coerced to ``str`` before grouping.
    layer : str, optional
        Name of an ``adata.layers`` entry to read expression from. If ``None``
        (default) or the name is absent, ``adata.X`` (or ``adata.raw.X`` when
        ``use_raw=True``) is used instead.
    use_raw : bool, default ``False``
        If ``True`` and ``adata.raw`` is present, read expression from
        ``adata.raw`` instead of ``adata``.
    expression_cutoff : float, default ``0``
        Threshold used to compute the expressing fraction (dot size): a cell
        counts as expressing a gene when its value is strictly greater than
        this cutoff.
    cmap : str, default ``"Reds"``
        Matplotlib colormap name mapping mean expression to dot colour.
    palette : dict, str or None, optional
        Colour spec for the primary ``Group`` (``groupby``) annotation stripe.
        Accepts an explicit ``{category: colour}`` mapping, a path to an Excel
        palette sheet, or ``None`` to auto-resolve from
        ``adata.uns[f"{groupby}_colors"]`` (falling back to matplotlib
        ``tab10``/``tab20``). See :func:`resolve_palette`.
    annotate : sequence of str, optional
        Extra ``adata.obs`` columns drawn as additional left-side row
        annotation stripes (one stripe per column). For each group the dominant
        (modal) value of the column is used.
    annotate_palettes : dict, optional
        Per-annotation-column palette spec, mapping an ``annotate`` column name
        to a palette (dict / Excel path / ``None``) resolved the same way as
        ``palette``.
    row_cluster : bool, default ``False``
        Whether to hierarchically cluster and reorder the rows (groups).
    col_cluster : bool, default ``False``
        Whether to hierarchically cluster and reorder the columns (genes).
    figsize : tuple of float, optional
        Figure size in inches. If ``None`` (default) it is derived from the
        number of genes and groups so dots and tick labels stay readable on
        large panels.
    title : str, optional
        Figure suptitle drawn above the dot plot. No title if ``None``.
    save : str, optional
        Path to write the figure to (parent directories created as needed).
        If ``None`` (default) the figure is not saved.
    show : bool, default ``False``
        If ``True`` display the figure with ``plt.show()``.
    show_rownames : bool, default ``True``
        Whether to draw row (group) labels.
    show_colnames : bool, default ``True``
        Whether to draw column (gene) labels.
    legend_kws : dict, optional
        Extra keyword arguments forwarded as the ``cmap_legend_kws`` of the
        underlying ``DotClustermapPlotter`` (colour legend styling).
    plot_kws : dict, optional
        Extra keyword arguments forwarded to
        :class:`PyComplexHeatmap.DotClustermapPlotter`, overriding the
        wrapper's defaults for any key they set.
    annot_kws : dict, optional
        Extra keyword arguments forwarded to
        :class:`PyComplexHeatmap.HeatmapAnnotation` used to build the row
        annotations.
    marker : str, default ``"o"``
        Matplotlib marker used for the dots.
    **kwargs
        Additional keyword arguments merged into ``plot_kws`` and passed
        straight through to ``DotClustermapPlotter``, so any PyComplexHeatmap
        parameter not exposed explicitly here remains reachable.

    Returns
    -------
    PyComplexHeatmap.DotClustermapPlotter
        The fitted plotter object holding the drawn dot heatmap (its ``fig`` /
        axes can be used for further customisation).
    """
    import PyComplexHeatmap as pch

    mean_df, frac_df = _aggregate(
        adata,
        groupby,
        list(genes),
        layer=layer,
        use_raw=use_raw,
        expression_cutoff=expression_cutoff,
    )
    cats = list(mean_df.index)
    long = (
        mean_df.stack()
        .rename("value")
        .reset_index()
        .rename(columns={mean_df.index.name or "level_0": groupby, "level_1": "gene"})
    )
    long["value"] = long["value"].astype(float)
    fl = (
        frac_df.stack()
        .rename("size")
        .reset_index()
        .rename(columns={frac_df.index.name or "level_0": groupby, "level_1": "gene"})
    )
    long["size"] = fl["size"].astype(float).values

    row_anno, _ = _row_annotation(
        adata, cats, groupby, palette, annotate, annotate_palettes, annot_kws
    )

    # Scale with gene (column) and group (row) counts so dots and tick
    # labels stay readable instead of overlapping on large panels.
    if figsize is None:
        n_col, n_row = len(mean_df.columns), len(mean_df.index)
        figsize = (max(4.0, 0.30 * n_col + 2.5), max(3.0, 0.28 * n_row + 1.8))

    fig = plt.figure(figsize=figsize)
    plot_kws = dict(plot_kws or {})
    plot_kws.update(kwargs)
    # Right-side row labels (row_names_side="right") extend past the heatmap
    # into the legend column; widen the legend gap for long group names so
    # the colour legend does not overlap the y tick labels.
    _maxlen = max((len(str(c)) for c in cats), default=0)
    _auto_hpad = max(4, round(0.75 * _maxlen))
    base = dict(
        data=long,
        x="gene",
        y=groupby,
        value="value",
        c="value",
        s="size",
        cmap=cmap,
        marker=marker,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        left_annotation=row_anno,
        verbose=0,
        show_rownames=show_rownames,
        show_colnames=show_colnames,
        col_names_side="bottom",
        row_names_side="right",
        xticklabels_kws=dict(labelrotation=-45, labelsize=8, bottom=True),
        yticklabels_kws=dict(labelsize=8, right=True),
        legend_hpad=_auto_hpad,
        legend_side="right",
        legend_anchor="ax_heatmap",
        spines=False,
    )
    for k, v in plot_kws.items():
        base[k] = v
    if legend_kws:
        base["cmap_legend_kws"] = dict(legend_kws)
    with _quiet_pch(fig):
        cm = pch.DotClustermapPlotter(**base)
    _strip_cbar_white_lines(cm)
    if title:
        fig.suptitle(title, y=1.02)
    with _quiet_pch(fig):
        save_or_show(fig, save, show=show)
    return cm
