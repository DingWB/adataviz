"""Composition-style plots: rose, ring, pie, area, dot, trend.

These visualise how categorical cell annotations (e.g. ``Class``,
``Subclass``) distribute across another categorical axis (e.g. donor,
brain region). All functions accept either an :class:`anndata.AnnData`
or a :class:`pandas.DataFrame`; when given an AnnData they automatically
pick up colour palettes from ``adata.uns[f"{groupby}_colors"]``.

Inspired by the ``CellStatPlot`` family of plot types, re-implemented in
pure matplotlib.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._utils import (
    categorical_order,
    maybe_adjust_texts,
    merge_legend_kws,
    resolve_adata_obs,
    resolve_palette,
    save_or_show,
)

__all__ = [
    "rose_plot",
    "ring_plot",
    "pie_plot",
    "area_plot",
    "dot_plot",
    "trend_plot",
    "bar_plot",
]


# ---------------------------------------------------------------------------
# Rose plot
# ---------------------------------------------------------------------------


def rose_plot(
    adata,
    groupby: str,
    split_by: Optional[str] = None,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=(5, 5),
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    edgecolor: str = "white",
    linewidth: float = 0.5,
    show_legend: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
    label_rotation: bool = True,
):
    """Draw a Nightingale rose (polar) chart of category counts.

    Parameters
    ----------
    adata : AnnData, DataFrame, or path
        Cell metadata source. Operates on ``adata.obs`` when given an AnnData.
    groupby : str
        Column whose categories define wedges.
    split_by : str, optional
        Column whose categories stack within each wedge.
    legend_kws : dict, optional
        Forwarded to :meth:`matplotlib.axes.Axes.legend`. Defaults follow
        ``plot_genes`` (right of axes, frameless, fontsize 8).
    label_rotation : bool, default True
        Rotate wedge labels to follow the radial direction (avoids
        overlap when there are many categories).
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    n = len(cats)
    angles = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    width = 2 * np.pi / max(n, 1) * 0.95

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection="polar")
    else:
        fig = ax.figure

    if split_by is None:
        colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
        ax.bar(
            angles,
            counts.values,
            width=width,
            bottom=0.0,
            color=[colors[c] for c in cats],
            edgecolor=edgecolor,
            linewidth=linewidth,
        )
    else:
        stack_cats = categorical_order(obs[split_by], split_order)
        ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
        ct = ct.reindex(index=cats, columns=stack_cats, fill_value=0)
        colors = resolve_palette(
            palette, stack_cats, sheet_name=split_by, adata=ad, groupby=split_by
        )
        bottom = np.zeros(n, dtype=float)
        for sc in stack_cats:
            vals = ct[sc].values.astype(float)
            ax.bar(
                angles,
                vals,
                width=width,
                bottom=bottom,
                color=colors[sc],
                edgecolor=edgecolor,
                linewidth=linewidth,
                label=sc,
            )
            bottom += vals
        if show_legend:
            ax.legend(**merge_legend_kws(legend_kws, title=split_by))

    ax.set_xticks(angles)
    if label_rotation and n > 6:
        ax.set_xticklabels([])
        for ang, lab in zip(angles, cats):
            rot = np.rad2deg(ang) - 90 if np.cos(ang) >= 0 else np.rad2deg(ang) + 90
            ax.text(
                ang,
                ax.get_rmax() * 1.05,
                lab,
                rotation=rot,
                ha="center",
                va="center",
                fontsize=8,
            )
    else:
        ax.set_xticklabels(cats, fontsize=8)
    ax.set_yticklabels([])
    ax.spines["polar"].set_visible(False)
    if title:
        ax.set_title(title)
    save_or_show(fig, save)
    return ax


# ---------------------------------------------------------------------------
# Ring plot
# ---------------------------------------------------------------------------


def ring_plot(
    adata,
    groupby: str,
    order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    width: float = 0.35,
    figsize=(4.5, 4.5),
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    autopct: Optional[str] = "%.1f%%",
    edgecolor: str = "white",
    show_legend: bool = False,
    legend_kws: Optional[Mapping[str, Any]] = None,
    show_labels: Optional[bool] = None,
    adjust_text: bool = False,
):
    """Donut chart for ``groupby`` value counts.

    With many categories (>= 12) labels around the ring overlap; pass
    ``show_labels=False`` (auto when there are many slices) to drop them
    and rely on ``show_legend=True`` instead, or pass ``adjust_text=True``
    to nudge them apart with :mod:`adjustText`.
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)

    if show_labels is None:
        show_labels = len(cats) <= 12

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    wedges, texts, autotexts = ax.pie(
        counts.values,
        labels=cats if show_labels else None,
        colors=[colors[c] for c in cats],
        wedgeprops=dict(width=width, edgecolor=edgecolor),
        autopct=autopct,
        pctdistance=1 - width / 2,
        textprops=dict(fontsize=8),
    )
    if adjust_text and show_labels:
        maybe_adjust_texts(texts, ax=ax)
    ax.set(aspect="equal")
    if show_legend:
        ax.legend(
            wedges,
            cats,
            **merge_legend_kws(legend_kws, title=groupby),
        )
    if title:
        ax.set_title(title)
    save_or_show(fig, save)
    return ax


# ---------------------------------------------------------------------------
# Pie plot
# ---------------------------------------------------------------------------


def pie_plot(
    adata,
    groupby: str,
    split_by: Optional[str] = None,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=None,
    ncols: Optional[int] = None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    autopct: Optional[str] = "%.1f%%",
    show_legend: Optional[bool] = None,
    legend_kws: Optional[Mapping[str, Any]] = None,
    show_labels: Optional[bool] = None,
    adjust_text: bool = False,
):
    """Pie chart of ``groupby`` counts, optionally faceted by ``split_by``.

    For a single pie returns the :class:`Axes`; for a faceted grid
    returns the array of axes. With many categories a shared legend is
    drawn and individual labels are suppressed automatically (override
    via ``show_labels``).
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
    if show_labels is None:
        show_labels = len(cats) <= 10
    if show_legend is None:
        show_legend = not show_labels

    if split_by is None:
        if figsize is None:
            figsize = (4.5, 4.5)
        fig, ax = plt.subplots(figsize=figsize)
        counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
        wedges, texts, _autotexts = ax.pie(
            counts.values,
            labels=cats if show_labels else None,
            colors=[colors[c] for c in cats],
            autopct=autopct,
            textprops=dict(fontsize=8),
        )
        if adjust_text and show_labels:
            maybe_adjust_texts(texts, ax=ax)
        if show_legend:
            ax.legend(wedges, cats, **merge_legend_kws(legend_kws, title=groupby))
        ax.set(aspect="equal")
        if title:
            ax.set_title(title)
        save_or_show(fig, save)
        return ax

    splits = categorical_order(obs[split_by], split_order)
    n = len(splits)
    if ncols is None:
        ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    if figsize is None:
        figsize = (3.2 * ncols, 3.2 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_1d(axes).ravel()
    for i, sp in enumerate(splits):
        sub = obs[obs[split_by].astype(str) == sp]
        counts = sub[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
        ax = axes[i]
        wedges, texts, _ap = ax.pie(
            counts.values,
            labels=cats if show_labels else None,
            colors=[colors[c] for c in cats],
            autopct=autopct,
            textprops=dict(fontsize=7),
        )
        if adjust_text and show_labels:
            maybe_adjust_texts(texts, ax=ax)
        ax.set_title(str(sp), fontsize=9)
        ax.set(aspect="equal")
    for j in range(n, len(axes)):
        axes[j].axis("off")
    if show_legend:
        from matplotlib.patches import Patch

        handles = [Patch(color=colors[c], label=c) for c in cats]
        fig.legend(
            handles=handles,
            **merge_legend_kws(
                legend_kws,
                loc="center left",
                bbox_to_anchor=(1.0, 0.5),
                title=groupby,
            ),
        )
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    save_or_show(fig, save)
    return axes


# ---------------------------------------------------------------------------
# Area plot
# ---------------------------------------------------------------------------


def area_plot(
    adata,
    groupby: str,
    split_by: str,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    normalize: bool = True,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=(6, 4),
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    alpha: float = 0.85,
    show_legend: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
):
    """Stacked area chart of ``groupby`` composition across ``split_by``."""
    obs, ad = resolve_adata_obs(adata)
    splits = categorical_order(obs[split_by], split_order)
    cats = categorical_order(obs[groupby], order)
    ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
    ct = ct.reindex(index=cats, columns=splits, fill_value=0).astype(float)
    if normalize:
        col_sums = ct.sum(axis=0).replace(0, np.nan)
        ct = ct.div(col_sums, axis=1).fillna(0)

    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    x = np.arange(len(splits))
    ax.stackplot(
        x,
        ct.values,
        labels=cats,
        colors=[colors[c] for c in cats],
        alpha=alpha,
        edgecolor="white",
        linewidth=0.3,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(splits, rotation=30, ha="right")
    ax.set_xlim(0, max(len(splits) - 1, 1))
    ax.set_ylim(0, 1 if normalize else None)
    ax.set_ylabel("Fraction" if normalize else "Count")
    ax.set_xlabel(split_by)
    if show_legend:
        ax.legend(**merge_legend_kws(legend_kws, title=groupby))
    if title:
        ax.set_title(title)
    save_or_show(fig, save)
    return ax


# ---------------------------------------------------------------------------
# Dot plot (count × column-fraction)
# ---------------------------------------------------------------------------


def dot_plot(
    adata,
    groupby: str,
    split_by: str,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    cmap: str = "viridis",
    size_max: float = 250.0,
    figsize=None,
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    cbar_label: Optional[str] = None,
):
    """Dot plot where size = count and colour = column-fraction.

    Rows are ``groupby`` categories, columns are ``split_by``. Useful
    for visualising the relative composition of one categorical variable
    across another.
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    splits = categorical_order(obs[split_by], split_order)
    ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
    ct = ct.reindex(index=cats, columns=splits, fill_value=0).astype(float)
    col_sums = ct.sum(axis=0).replace(0, np.nan)
    frac = ct.div(col_sums, axis=1).fillna(0).values
    counts = ct.values

    if figsize is None:
        figsize = (max(4, 0.5 * len(splits) + 2), max(3, 0.35 * len(cats) + 1.5))
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    yy, xx = np.meshgrid(np.arange(len(cats)), np.arange(len(splits)), indexing="ij")
    max_count = max(counts.max(), 1)
    sizes = (counts / max_count) * size_max
    sc = ax.scatter(
        xx.ravel(),
        yy.ravel(),
        s=sizes.ravel(),
        c=frac.ravel(),
        cmap=cmap,
        edgecolor="black",
        linewidth=0.3,
        vmin=0,
        vmax=1,
    )
    ax.set_xticks(np.arange(len(splits)))
    ax.set_xticklabels(splits, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(cats)))
    ax.set_yticklabels(cats)
    ax.set_xlim(-0.5, len(splits) - 0.5)
    ax.set_ylim(-0.5, len(cats) - 0.5)
    ax.invert_yaxis()
    ax.set_xlabel(split_by)
    ax.set_ylabel(groupby)
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.6)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.04)
    cbar.set_label(cbar_label or f"Fraction within {split_by}")
    if title:
        ax.set_title(title)
    save_or_show(fig, save)
    return ax


# ---------------------------------------------------------------------------
# Trend plot
# ---------------------------------------------------------------------------


def trend_plot(
    adata,
    groupby: str,
    split_by: str,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    normalize: bool = True,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=(6, 4),
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    marker: str = "o",
    show_legend: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
):
    """Line trend of each ``groupby`` category across ordered ``split_by``."""
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    splits = categorical_order(obs[split_by], split_order)
    ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
    ct = ct.reindex(index=cats, columns=splits, fill_value=0).astype(float)
    if normalize:
        col_sums = ct.sum(axis=0).replace(0, np.nan)
        ct = ct.div(col_sums, axis=1).fillna(0)

    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    x = np.arange(len(splits))
    for c in cats:
        ax.plot(
            x,
            ct.loc[c].values,
            marker=marker,
            color=colors[c],
            label=c,
            linewidth=1.5,
            markersize=4,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(splits, rotation=30, ha="right")
    ax.set_xlabel(split_by)
    ax.set_ylabel("Fraction" if normalize else "Count")
    if show_legend:
        ax.legend(**merge_legend_kws(legend_kws, title=groupby))
    if title:
        ax.set_title(title)
    save_or_show(fig, save)
    return ax


# ---------------------------------------------------------------------------
# Stacked bar plot
# ---------------------------------------------------------------------------


def bar_plot(
    adata,
    groupby: str,
    split_by: str,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    normalize: bool = True,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=None,
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    rotation: float = -45,
    width: float = 0.95,
    edgecolor: str = "black",
    linewidth: float = 0.1,
    show_legend: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
):
    """Stacked bar chart of ``groupby`` composition across ``split_by``.

    Each bar represents one ``split_by`` category and is stacked into
    ``groupby`` proportions (or absolute counts when ``normalize=False``).
    """
    obs, ad = resolve_adata_obs(adata)
    splits = categorical_order(obs[split_by], split_order)
    cats = categorical_order(obs[groupby], order)
    ct = pd.crosstab(obs[split_by].astype(str), obs[groupby].astype(str))
    ct = ct.reindex(index=splits, columns=cats, fill_value=0).astype(float)
    if normalize:
        row_sums = ct.sum(axis=1).replace(0, np.nan)
        ct = ct.div(row_sums, axis=0).fillna(0)
    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
    if figsize is None:
        figsize = (max(3, len(splits) * 0.45), 4.5)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    bottom = np.zeros(len(splits), dtype=float)
    x = np.arange(len(splits))
    for c in cats:
        vals = ct[c].values
        ax.bar(
            x, vals, width=width, bottom=bottom,
            color=colors[c], edgecolor=edgecolor, linewidth=linewidth,
            label=c,
        )
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(splits, rotation=rotation, ha="left" if rotation < 0 else "right")
    ax.set_xlabel(split_by)
    ax.set_ylabel("Fraction" if normalize else "Count")
    if normalize:
        ax.set_ylim(0, 1)
    if show_legend:
        ax.legend(**merge_legend_kws(legend_kws, title=groupby))
    if title:
        ax.set_title(title)
    save_or_show(fig, save, show=show)
    return ax
