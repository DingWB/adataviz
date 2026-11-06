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

from typing import Any, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._utils import (
    categorical_series,
    gene_positions,
    import_pch,
    quiet_pch as _quiet_pch,
    resolve_palette,
    save_or_show,
)

__all__ = ["complex_heatmap", "complex_dotplot"]

# PyComplexHeatmap's dot-size legend uses ``handlelength=1`` (one em), but the
# reference dots are as large as a heatmap cell (~20 pt) and overflow onto
# their labels; give the handles room and space the rows out.
_DOT_LEGEND_KWS = dict(handlelength=3.0, handletextpad=0.8, labelspacing=0.9)


def _aggregate(adata, groupby, genes, layer=None, use_raw=False, expression_cutoff=0):
    """Return ``(mean_df, frac_df)`` of ``genes`` per ``groupby`` category.

    Groups keep the categorical order of ``adata.obs[groupby]`` (or first
    appearance order otherwise) so that un-clustered heatmaps follow the
    user's category ordering instead of being sorted alphabetically.
    """
    import anndata as _ad

    if not isinstance(adata, _ad.AnnData):
        raise TypeError("complex_heatmap/complex_dotplot require AnnData input.")
    src = adata.raw if use_raw and adata.raw is not None else adata
    # Column positions must come from the matrix we are about to slice: with
    # use_raw=True, adata.raw.var_names and adata.var_names can differ in both
    # order and length, so indices built from one silently select the wrong
    # genes out of the other.
    if layer is not None and layer in adata.layers:
        mat = adata.layers[layer]
        var_names = adata.var_names
    else:
        mat = src.X
        var_names = src.var_names
    X = mat[:, gene_positions(var_names, genes)]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)
    g = categorical_series(adata.obs[groupby])
    df = pd.DataFrame(X, index=adata.obs_names, columns=list(genes))
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
    pch = import_pch()

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
        g = categorical_series(adata.obs[groupby])
        for col in annotate:
            if col not in adata.obs.columns:
                continue
            # Dominant (modal) value of ``col`` per group, without a python
            # lambda per group: count (group, value) pairs once and keep the
            # most frequent value for each group.
            pairs = pd.DataFrame({"g": g, "v": adata.obs[col].astype(str).values})
            modal = (
                pairs.groupby(["g", "v"], observed=True)
                .size()
                .reset_index(name="n")
                .sort_values(["g", "n"], ascending=[True, False])
                .drop_duplicates("g")
                .set_index("g")["v"]
                .reindex(cats)
            )
            ann_cats = modal.astype(str).tolist()
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


def _auto_figsize(n_col: int, n_row: int, cell: float = 0.28, extra_w: float = 2.5):
    """Figure size that grows with the grid so cells stay roughly square."""
    return (max(4.0, cell * n_col + extra_w), max(3.0, cell * n_row + 1.8))


def _set_title(fig, cm, title: Optional[str]):
    """Put ``title`` centred above the heatmap/annotation block.

    ``fig.suptitle`` centres on the whole figure (which includes the legend
    column) so it ends up off-centre; align it with the heatmap axes instead
    and sit it just above the topmost axes of the plotter.
    """
    if not title:
        return
    try:
        fig.canvas.draw()
        hm = cm.ax_heatmap.get_position()
        top = max(a.get_position().y1 for a in fig.axes if a.get_visible())
        fig.text(
            (hm.x0 + hm.x1) / 2,
            top + 0.02,
            title,
            ha="center",
            va="bottom",
            fontsize=plt.rcParams["figure.titlesize"],
        )
    except Exception:
        fig.suptitle(title, y=1.02)


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
    pch = import_pch()

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
        figsize = _auto_figsize(len(mean_df.columns), len(mean_df.index))

    fig = plt.figure(figsize=figsize)
    plot_kws = dict(plot_kws or {})
    plot_kws.update(kwargs)
    base = dict(
        data=mean_df,
        cmap=cmap,
        # Colorbar title: say what the colours mean instead of PCH's "heatmap".
        label={"row": "z-score (row)", "col": "z-score (col)"}.get(
            z_score, "Mean expression"
        ),
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        show_rownames=show_rownames,
        show_colnames=show_colnames,
        left_annotation=row_anno,
        col_names_side="bottom",
        row_names_side="right",
        xticklabels_kws=dict(labelrotation=-45, labelsize=8, bottom=True),
        yticklabels_kws=dict(labelsize=8, right=True),
        # PCH already reserves the row-label width when row names sit on the
        # right; only a small extra gap is needed before the legends.
        legend_hpad=3,
        legend_side="right",
        legend_anchor="ax_heatmap",
        verbose=0,
    )
    base.update(plot_kws)
    if legend_kws:
        base["legend_kws"] = dict(legend_kws)
    with _quiet_pch(fig):
        cm = pch.ClusterMapPlotter(**base)
    _strip_cbar_white_lines(cm)
    _set_title(fig, cm, title)
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
    pch = import_pch()

    mean_df, frac_df = _aggregate(
        adata,
        groupby,
        list(genes),
        layer=layer,
        use_raw=use_raw,
        expression_cutoff=expression_cutoff,
    )
    cats = list(mean_df.index)
    gene_list = list(mean_df.columns)
    # Long table straight from the aligned matrices (same index/columns), so
    # no reliance on two independent ``stack()`` calls lining up row-by-row.
    long = pd.DataFrame(
        {
            groupby: np.repeat(np.asarray(cats, dtype=object), len(gene_list)),
            "gene": np.tile(np.asarray(gene_list, dtype=object), len(cats)),
            "Mean": mean_df.to_numpy(dtype=float).ravel(),
            "Fraction": frac_df.to_numpy(dtype=float).ravel(),
        }
    )

    row_anno, _ = _row_annotation(
        adata, cats, groupby, palette, annotate, annotate_palettes, annot_kws
    )

    # Scale with gene (column) and group (row) counts so dots and tick
    # labels stay readable instead of overlapping on large panels.
    if figsize is None:
        figsize = _auto_figsize(len(gene_list), len(cats), cell=0.30)

    fig = plt.figure(figsize=figsize)
    plot_kws = dict(plot_kws or {})
    plot_kws.update(kwargs)
    cmap_legend_kws = dict(label="Mean expression", extend="both", extendfrac=0.1)
    if legend_kws:
        cmap_legend_kws.update(legend_kws)
    base = dict(
        data=long,
        x="gene",
        y=groupby,
        value="Mean",
        c="Mean",
        s="Fraction",
        cmap=cmap,
        marker=marker,
        dot_legend_marker=marker,
        cmap_legend_kws=cmap_legend_kws,
        dot_legend_kws=dict(_DOT_LEGEND_KWS),
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        # Honour the categorical / requested order when not clustering.
        x_order=None if col_cluster else gene_list,
        y_order=None if row_cluster else cats,
        left_annotation=row_anno,
        verbose=0,
        show_rownames=show_rownames,
        show_colnames=show_colnames,
        col_names_side="bottom",
        row_names_side="right",
        xticklabels_kws=dict(labelrotation=-45, labelsize=8, bottom=True),
        yticklabels_kws=dict(labelsize=8, right=True),
        legend_hpad=3,
        legend_side="right",
        legend_anchor="ax_heatmap",
        spines=False,
    )
    base.update(plot_kws)
    with _quiet_pch(fig):
        cm = pch.DotClustermapPlotter(**base)
    _strip_cbar_white_lines(cm)
    _set_title(fig, cm, title)
    with _quiet_pch(fig):
        save_or_show(fig, save, show=show)
    return cm
