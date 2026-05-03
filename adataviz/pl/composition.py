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
    label_orientation: str = "radial",
    label_pad: float = 1.12,
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
        Forwarded to :meth:`matplotlib.axes.Axes.legend`.
    label_orientation : {"radial", "vertical", "none"}, default ``"radial"``
        ``"radial"`` rotates each label tangent to its wedge (flipped on the
        lower half so it stays right-side-up); ``"vertical"`` keeps labels
        upright; ``"none"`` hides labels.
    label_pad : float, default 1.12
        Multiplier on ``rmax`` controlling label distance from the outer
        ring (raise to avoid title overlap).
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
    if label_orientation == "none":
        ax.set_xticklabels([])
    elif label_orientation == "radial":
        ax.set_xticklabels([])
        rmax = ax.get_rmax()
        for ang, lab in zip(angles, cats):
            deg = (np.rad2deg(ang) + 360) % 360
            if deg <= 90 or deg >= 270:
                rot, ha = deg, "left"
            else:
                rot, ha = deg - 180, "right"
            ax.text(
                ang, rmax * label_pad, lab,
                rotation=rot, rotation_mode="anchor",
                ha=ha, va="center", fontsize=8,
            )
    else:  # vertical / upright
        ax.set_xticklabels(cats, fontsize=8)
    ax.set_yticklabels([])
    ax.spines["polar"].set_visible(False)
    ax.grid(color="lightgrey", linewidth=0.4, alpha=0.7)
    # Reserve enough margin so labels don't collide with the title.
    if title:
        ax.set_title(title, pad=24)
    fig.subplots_adjust(top=0.86, left=0.05, right=0.95, bottom=0.05)
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
    legend_loc: str = "right",
    show_labels: Optional[bool] = None,
    adjust_text: bool = False,
    label_orientation: str = "horizontal",
):
    """Donut chart for ``groupby`` value counts.

    Parameters
    ----------
    legend_loc : {"right", "center"}, default ``"right"``
        ``"right"`` places the legend to the right of the donut
        (useful when there are many categories); ``"center"`` draws it
        inside the donut hole - clean look when categories are few.
    label_orientation : {"horizontal", "radial"}, default ``"horizontal"``
        ``"radial"`` rotates each external label tangent to its wedge
        (flipped on the lower half so it stays readable).
    show_labels : bool, optional
        Auto-disabled when there are >12 categories. Pass
        ``adjust_text=True`` to nudge labels apart with :mod:`adjustText`.
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)

    if show_labels is None:
        show_labels = len(cats) <= 12 and legend_loc != "center"

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    total = float(counts.sum()) or 1.0
    fracs = counts.values / total
    # Hide percent labels on tiny slices to avoid overlap.
    def _autopct(p):
        return (autopct % p) if autopct and p >= 3 else ""
    wedges, texts, autotexts = ax.pie(
        counts.values,
        labels=cats if show_labels else None,
        colors=[colors[c] for c in cats],
        wedgeprops=dict(width=width, edgecolor=edgecolor, linewidth=0.6),
        autopct=_autopct if autopct else None,
        pctdistance=1 - width / 2,
        textprops=dict(fontsize=7),
    )
    # Drop labels for tiny slices to declutter.
    if show_labels:
        for t, f in zip(texts, fracs):
            if f < 0.02:
                t.set_text("")
        if label_orientation == "radial":
            import numpy as _np

            for t, w in zip(texts, wedges):
                if not t.get_text():
                    continue
                ang = (w.theta1 + w.theta2) / 2.0
                deg = (ang + 360) % 360
                if deg <= 90 or deg >= 270:
                    rot, ha = deg, "left"
                else:
                    rot, ha = deg - 180, "right"
                # Position just outside the donut.
                r = 1.05
                t.set_position((r * _np.cos(_np.deg2rad(ang)), r * _np.sin(_np.deg2rad(ang))))
                t.set_rotation(rot)
                t.set_rotation_mode("anchor")
                t.set_ha(ha)
                t.set_va("center")
        if adjust_text:
            maybe_adjust_texts([t for t in texts if t.get_text()], ax=ax)
    ax.set(aspect="equal")
    if show_legend:
        if legend_loc == "center":
            lkws = dict(
                loc="center", bbox_to_anchor=(0.5, 0.5), frameon=False,
                fontsize=8, title=groupby, title_fontsize=9,
                handlelength=1.0, handletextpad=0.4, labelspacing=0.3,
            )
            if legend_kws:
                lkws.update(legend_kws)
            ax.legend(wedges, cats, **lkws)
        else:
            ax.legend(wedges, cats, **merge_legend_kws(legend_kws, title=groupby))
    if title:
        ax.set_title(title, pad=10)
    fig.tight_layout()
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
        total = float(counts.sum()) or 1.0
        fracs = counts.values / total
        def _autopct(p):
            return (autopct % p) if autopct and p >= 3 else ""
        wedges, texts, _autotexts = ax.pie(
            counts.values,
            labels=cats if show_labels else None,
            colors=[colors[c] for c in cats],
            autopct=_autopct if autopct else None,
            textprops=dict(fontsize=8),
            wedgeprops=dict(edgecolor="white", linewidth=0.6),
        )
        if show_labels:
            for t, f in zip(texts, fracs):
                if f < 0.02:
                    t.set_text("")
            if adjust_text:
                maybe_adjust_texts([t for t in texts if t.get_text()], ax=ax)
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
        total = float(counts.sum()) or 1.0
        fracs = counts.values / total
        def _autopct_facet(p):
            return (autopct % p) if autopct and p >= 3 else ""
        wedges, texts, _ap = ax.pie(
            counts.values,
            labels=cats if show_labels else None,
            colors=[colors[c] for c in cats],
            autopct=_autopct_facet if autopct else None,
            textprops=dict(fontsize=7),
            wedgeprops=dict(edgecolor="white", linewidth=0.5),
        )
        if show_labels:
            for t, f in zip(texts, fracs):
                if f < 0.02:
                    t.set_text("")
            if adjust_text:
                maybe_adjust_texts([t for t in texts if t.get_text()], ax=ax)
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
    ax.set_xticklabels(splits, rotation=-45, rotation_mode="anchor", ha="left", va="center")
    ax.set_xlim(0, max(len(splits) - 1, 1))
    ax.set_ylim(0, 1 if normalize else None)
    ax.set_ylabel("Fraction" if normalize else "Count")
    ax.set_xlabel(split_by)
    ax.spines[["top", "right"]].set_visible(False)
    if show_legend:
        ax.legend(**merge_legend_kws(legend_kws, title=groupby))
    if title:
        ax.set_title(title)
    fig.tight_layout()
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
    ax.set_xticklabels(splits, rotation=-45, rotation_mode="anchor", ha="left", va="center")
    ax.set_yticks(np.arange(len(cats)))
    ax.set_yticklabels(cats)
    ax.set_xlim(-0.5, len(splits) - 0.5)
    ax.set_ylim(-0.5, len(cats) - 0.5)
    ax.invert_yaxis()
    ax.set_xlabel(split_by)
    ax.set_ylabel(groupby)
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.04)
    cbar.set_label(cbar_label or f"Fraction within {split_by}", fontsize=8)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=2, width=0.5, labelsize=7)
    if title:
        ax.set_title(title)
    fig.tight_layout()
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
    ax.set_xticklabels(splits, rotation=-45, rotation_mode="anchor", ha="left", va="center")
    ax.set_xlim(-0.3, len(splits) - 0.7)
    ax.set_xlabel(split_by)
    ax.set_ylabel("Fraction" if normalize else "Count")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linewidth=0.4, alpha=0.5)
    if show_legend:
        ax.legend(**merge_legend_kws(legend_kws, title=groupby))
    if title:
        ax.set_title(title)
    fig.tight_layout()
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
    gap: float = 0.05,
    edgecolor: str = "black",
    linewidth: float = 0.1,
    show_legend: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
    sort_by: Optional[str] = None,
):
    """Stacked bar chart of ``groupby`` composition across ``split_by``.

    Each bar represents one ``split_by`` category and is stacked into
    ``groupby`` proportions (or absolute counts when ``normalize=False``).
    The y-axis is intentionally hidden (no ticks / labels / spine) for a
    cleaner publication-style composition plot, mirroring the legacy
    ``stacked_barplot`` API.
    """
    obs, ad = resolve_adata_obs(adata)
    splits = categorical_order(obs[split_by], split_order)
    cats = categorical_order(obs[groupby], order)
    ct = pd.crosstab(obs[split_by].astype(str), obs[groupby].astype(str))
    ct = ct.reindex(index=splits, columns=cats, fill_value=0).astype(float)
    if normalize:
        row_sums = ct.sum(axis=1).replace(0, np.nan)
        ct = ct.div(row_sums, axis=0).fillna(0)
    if sort_by is not None and sort_by in ct.columns:
        ct = ct.sort_values(sort_by, ascending=True)
        splits = list(ct.index)
    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
    if figsize is None:
        figsize = (max(2.5, len(splits) * 0.45), max(3.5, min(len(cats) * 0.5, 8)))
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    n = len(splits)
    bottom = np.zeros(n, dtype=float)
    x = np.arange(n, dtype=float)  # left edges
    bar_width = 1.0 - gap
    for c in cats:
        vals = ct[c].values
        ax.bar(
            x, vals, width=bar_width, bottom=bottom,
            align="edge",
            color=colors[c], edgecolor=edgecolor, linewidth=linewidth,
            label=c,
        )
        bottom += vals
    # Tight x-limits remove the empty matplotlib auto-margin on either side.
    ax.set_xlim(0, n)
    if normalize:
        ax.set_ylim(0, 1)
    ax.set_xticks(x + 0.5)
    ax.set_xticklabels(
        splits,
        rotation=rotation,
        rotation_mode="anchor",
        ha="left" if rotation < 0 else "right",
        va="center",
    )
    ax.tick_params(
        axis="y", which="both",
        left=False, right=False,
        labelleft=False, labelright=False,
    )
    ax.tick_params(axis="x", length=3, pad=2, top=False, labeltop=False)
    ax.xaxis.label.set_visible(False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(False)
    if show_legend:
        lkws = dict(
            loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False,
            fontsize=8, ncol=1,
            borderpad=0.3, handlelength=1.0, handleheight=1.0,
            handletextpad=0.4, labelspacing=0.25, columnspacing=0.8,
            title=groupby, title_fontsize=9,
        )
        if legend_kws:
            lkws.update(legend_kws)
        ax.legend(**lkws)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    save_or_show(fig, save, show=show)
    return ax
