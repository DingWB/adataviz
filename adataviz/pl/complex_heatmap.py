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

from ._utils import resolve_palette, save_or_show

__all__ = ["complex_heatmap", "complex_dotplot"]


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

    Parameters
    ----------
    annotate : list of str, optional
        Extra ``adata.obs`` columns drawn as additional row annotations
        (one stripe per column). The dominant value per ``groupby`` is used.
    annotate_palettes : dict, optional
        Per-annotation-column palette spec.
    plot_kws : dict, optional
        Forwarded to :class:`PyComplexHeatmap.ClusterMapPlotter`.
    annot_kws : dict, optional
        Forwarded to :class:`PyComplexHeatmap.HeatmapAnnotation`.
    **kwargs
        Any extra keyword goes straight into ``ClusterMapPlotter``.
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
        legend_hpad=4,
        legend_side="right",
        legend_anchor="ax_heatmap",
        verbose=0,
    )
    for k, v in plot_kws.items():
        base[k] = v
    if legend_kws:
        base["legend_kws"] = dict(legend_kws)
    cm = pch.ClusterMapPlotter(**base)
    _strip_cbar_white_lines(cm)
    if title:
        fig.suptitle(title, y=1.02)
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
    """Dot heatmap: dot size = fraction expressing, colour = mean expression."""
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
        legend_hpad=4,
        legend_side="right",
        legend_anchor="ax_heatmap",
        spines=False,
    )
    for k, v in plot_kws.items():
        base[k] = v
    if legend_kws:
        base["cmap_legend_kws"] = dict(legend_kws)
    cm = pch.DotClustermapPlotter(**base)
    _strip_cbar_white_lines(cm)
    if title:
        fig.suptitle(title, y=1.02)
    save_or_show(fig, save, show=show)
    return cm
