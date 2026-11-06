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
    apply_layout,
    categorical_order,
    label_extent,
    maybe_adjust_texts,
    outside_legend,
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


def _legend_aware_figsize(base, cats, show_legend: bool, fontsize: float = 8):
    """Grow ``base`` (w, h) so an outside legend for ``cats`` fits next to the axes.

    The legend column needs roughly the longest label width plus a handle,
    and the figure must be at least as tall as the stacked legend entries;
    otherwise matplotlib squeezes the axes to a sliver (or lets the legend
    collide with the tick labels).
    """
    w, h = base
    if not show_legend:
        return (w, h)
    w += label_extent(cats, fontsize) + 0.4
    h = max(h, 0.19 * len(cats) + 0.8)
    return (w, h)


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

    Each ``groupby`` category is drawn as a wedge on a polar axis whose
    radial length encodes its cell count. When ``split_by`` is given, every
    wedge is stacked into the ``split_by`` sub-categories, giving a compact
    circular view of how one categorical annotation is composed across
    another.

    Parameters
    ----------
    adata : AnnData, pandas.DataFrame, or str
        Cell metadata source. When an :class:`~anndata.AnnData` is passed the
        function operates on ``adata.obs`` and can pick up colours from
        ``adata.uns[f"{groupby}_colors"]``; a DataFrame or a path to one is
        also accepted.
    groupby : str
        Column in the metadata whose categories define the wedges (one wedge
        per category, sized by its count).
    split_by : str, optional
        Column whose categories are stacked within each ``groupby`` wedge. If
        ``None`` (default), each wedge is a single solid bar.
    order : sequence, optional
        Explicit ordering of the ``groupby`` categories around the circle.
        Defaults to the natural categorical order of the column.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` stack segments. Only used when
        ``split_by`` is provided.
    palette : dict, str, or None, default None
        Colour mapping. A dict maps category -> colour; a str is treated as a
        path to a palette workbook/sheet; ``None`` (default) resolves colours
        from the AnnData ``uns`` colours or a generated default. When
        ``split_by`` is given the palette applies to the ``split_by``
        categories, otherwise to the ``groupby`` categories.
    figsize : tuple of float, default (5, 5)
        Figure size in inches, used only when a new figure is created
        (i.e. ``ax`` is ``None``).
    ax : matplotlib.axes.Axes, optional
        Existing polar axes to draw into. If ``None`` (default) a new figure
        with a polar subplot is created.
    save : str, optional
        If given, the figure is written to this path.
    show : bool, default False
        Whether to display the figure interactively.
    title : str, optional
        Axes title. Extra top padding is reserved so it clears the labels.
    edgecolor : str, default "white"
        Colour of the wedge borders.
    linewidth : float, default 0.5
        Width of the wedge borders.
    show_legend : bool, default True
        Whether to draw a legend. Only relevant when ``split_by`` is set
        (a single-category rose has no legend).
    legend_kws : dict, optional
        Extra keyword arguments forwarded to
        :meth:`matplotlib.axes.Axes.legend` (merged with sensible defaults).
    label_orientation : {"radial", "vertical", "none"}, default "radial"
        ``"radial"`` rotates each label tangent to its wedge (flipped on the
        lower half so it stays right-side-up); ``"vertical"`` keeps labels
        upright; ``"none"`` hides labels.
    label_pad : float, default 1.12
        Multiplier on ``rmax`` controlling the label distance from the outer
        ring (raise to avoid title overlap).

    Returns
    -------
    matplotlib.axes.Axes
        The polar axes containing the rose plot.
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    n = len(cats)
    angles = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    width = 2 * np.pi / max(n, 1) * 0.95

    own_fig = ax is None
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection="polar")
    else:
        fig = ax.figure
    # Constrained layout shrinks the polar axes so that the category labels
    # placed outside the ring (and the legend / title) all fit the figure.
    apply_layout(fig, own=own_fig)

    if split_by is None:
        colors = resolve_palette(
            palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
        )
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
            # Figure-level legend: an axes-anchored one would sit on top of
            # the labels that stick out past the polar axes on the right.
            outside_legend(fig, ax, legend_kws=legend_kws, title=split_by)

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
                ang,
                rmax * label_pad,
                lab,
                rotation=rot,
                rotation_mode="anchor",
                ha=ha,
                va="center",
                fontsize=8,
            )
    else:  # vertical / upright
        ax.set_xticklabels(cats, fontsize=8)
        # Push upright labels outward on the side they sit on so they do
        # not overlap the outer ring of wedges.
        for lab, ang in zip(ax.get_xticklabels(), angles):
            c = np.cos(ang)
            lab.set_ha("left" if c > 0.1 else "right" if c < -0.1 else "center")
    ax.set_yticklabels([])
    ax.spines["polar"].set_visible(False)
    ax.grid(color="lightgrey", linewidth=0.4, alpha=0.7)
    if title:
        # suptitle sits above the outward-pointing labels; an axes title
        # would be drawn right across the labels at the top of the ring.
        if own_fig:
            fig.suptitle(title)
        else:
            ax.set_title(title, pad=24)
    save_or_show(fig, save, show=show)
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
    figsize=None,
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
    """Donut (ring) chart of ``groupby`` value counts.

    Draws a single pie with a hollow centre so the proportions of each
    ``groupby`` category read as ring segments. Percentage labels on very
    small slices are suppressed to reduce clutter, and category labels can
    optionally be placed around the ring or replaced by a legend.

    Parameters
    ----------
    adata : AnnData, pandas.DataFrame, or str
        Cell metadata source. When an :class:`~anndata.AnnData` is passed the
        function operates on ``adata.obs`` and can pick up colours from
        ``adata.uns[f"{groupby}_colors"]``; a DataFrame or a path to one is
        also accepted.
    groupby : str
        Column whose category counts are turned into ring segments.
    order : sequence, optional
        Explicit ordering of the ``groupby`` categories. Defaults to the
        natural categorical order of the column.
    palette : dict, str, or None, default None
        Colour mapping. A dict maps category -> colour; a str is treated as a
        path to a palette workbook/sheet; ``None`` (default) resolves colours
        from the AnnData ``uns`` colours or a generated default.
    width : float, default 0.35
        Radial width of the ring as a fraction of the radius (smaller values
        yield a thinner donut with a larger hole).
    figsize : tuple of float, optional
        Figure size in inches, used only when a new figure is created
        (``ax`` is ``None``). When ``None`` (default) a 4.5 inch square is
        widened / heightened as needed for an outside legend.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw into. If ``None`` (default) a new figure and
        axes are created.
    save : str, optional
        If given, the figure is written to this path.
    show : bool, default False
        Whether to display the figure interactively.
    title : str, optional
        Axes title.
    autopct : str or None, default "%.1f%%"
        Percentage format string for slice labels. Slices below 3%% are left
        unlabelled; pass ``None`` to disable percentage labels entirely.
    edgecolor : str, default "white"
        Colour of the wedge borders.
    show_legend : bool, default False
        Whether to draw a legend of the categories.
    legend_kws : dict, optional
        Extra keyword arguments forwarded to
        :meth:`matplotlib.axes.Axes.legend` (merged with defaults).
    legend_loc : {"right", "center"}, default "right"
        ``"right"`` places the legend to the right of the donut (useful when
        there are many categories); ``"center"`` draws it inside the donut
        hole - a clean look when categories are few.
    show_labels : bool, optional
        Whether to draw category labels around the ring. When ``None``
        (default) this is auto-enabled only if there are <=12 categories and
        ``legend_loc`` is not ``"center"``. Labels on slices below 2%% are
        dropped to declutter.
    adjust_text : bool, default False
        If ``True``, nudge overlapping labels apart with :mod:`adjustText`.
    label_orientation : {"horizontal", "radial"}, default "horizontal"
        ``"radial"`` rotates each external label tangent to its wedge
        (flipped on the lower half so it stays readable); ``"horizontal"``
        keeps labels upright.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the ring plot.
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )

    if show_labels is None:
        show_labels = len(cats) <= 12 and legend_loc != "center"
    if figsize is None:
        figsize = _legend_aware_figsize(
            (4.5, 4.5), cats, show_legend and legend_loc != "center"
        )

    own_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    apply_layout(fig, own=own_fig)

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
                t.set_position(
                    (r * np.cos(np.deg2rad(ang)), r * np.sin(np.deg2rad(ang)))
                )
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
                loc="center",
                bbox_to_anchor=(0.5, 0.5),
                frameon=False,
                fontsize=8,
                title=groupby,
                title_fontsize=9,
                handlelength=1.0,
                handletextpad=0.4,
                labelspacing=0.3,
            )
            if legend_kws:
                lkws.update(legend_kws)
            ax.legend(wedges, cats, **lkws)
        else:
            outside_legend(fig, ax, wedges, cats, legend_kws=legend_kws, title=groupby)
    if title:
        ax.set_title(title, pad=10)
    save_or_show(fig, save, show=show)
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

    Draws a standard pie of the ``groupby`` category proportions. When
    ``split_by`` is provided, a grid of pies is produced - one per
    ``split_by`` category - so the composition can be compared across
    conditions. With many categories a shared legend is drawn and per-slice
    labels are suppressed automatically (override via ``show_labels``).

    Parameters
    ----------
    adata : AnnData, pandas.DataFrame, or str
        Cell metadata source. When an :class:`~anndata.AnnData` is passed the
        function operates on ``adata.obs`` and can pick up colours from
        ``adata.uns[f"{groupby}_colors"]``; a DataFrame or a path to one is
        also accepted.
    groupby : str
        Column whose category counts define the pie slices.
    split_by : str, optional
        Column used to facet the plot into one pie per category. If ``None``
        (default) a single pie is drawn.
    order : sequence, optional
        Explicit ordering of the ``groupby`` slices. Defaults to the natural
        categorical order of the column.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` facets. Only used when
        ``split_by`` is provided.
    palette : dict, str, or None, default None
        Colour mapping. A dict maps category -> colour; a str is treated as a
        path to a palette workbook/sheet; ``None`` (default) resolves colours
        from the AnnData ``uns`` colours or a generated default.
    figsize : tuple of float, optional
        Figure size in inches. When ``None`` (default) a sensible size is
        chosen automatically: ``(4.5, 4.5)`` for a single pie (widened for
        an outside legend), or a size scaled by the facet grid dimensions.
    ncols : int, optional
        Number of columns in the facet grid. Only used when ``split_by`` is
        set; defaults to ``min(n_facets, 4)``.
    save : str, optional
        If given, the figure is written to this path.
    show : bool, default False
        Whether to display the figure interactively.
    title : str, optional
        Overall title. Used as the axes title for a single pie or the figure
        suptitle for a facet grid.
    autopct : str or None, default "%.1f%%"
        Percentage format string for slice labels. Slices below 4%% are left
        unlabelled; pass ``None`` to disable percentage labels entirely.
    show_legend : bool, optional
        Whether to draw a shared legend. When ``None`` (default) it is
        enabled automatically whenever per-slice labels are hidden
        (i.e. it is the inverse of ``show_labels``).
    legend_kws : dict, optional
        Extra keyword arguments forwarded to the legend (merged with
        defaults).
    show_labels : bool, optional
        Whether to draw slice labels. When ``None`` (default) labels are
        enabled only for a single pie with <=10 categories; facet grids
        (``split_by``) rely on the shared legend instead of repeating every
        label on every pie. Labels on slices below 2%% are dropped to
        declutter.
    adjust_text : bool, default False
        If ``True``, nudge overlapping labels apart with :mod:`adjustText`.

    Returns
    -------
    matplotlib.axes.Axes or numpy.ndarray of Axes
        A single axes when ``split_by`` is ``None``; otherwise the flattened
        array of facet axes.
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
    if show_labels is None:
        # Facet grids repeat every label on every pie: one shared legend is
        # much cleaner there, so labels default to off when ``split_by`` is set.
        show_labels = split_by is None and len(cats) <= 10
    if show_legend is None:
        show_legend = not show_labels

    def _autopct(p):
        # Hide percentages on thin slices so they can't pile onto each other.
        return (autopct % p) if autopct and p >= 4 else ""

    if split_by is None:
        if figsize is None:
            figsize = _legend_aware_figsize((4.5, 4.5), cats, show_legend)
        fig, ax = plt.subplots(figsize=figsize)
        apply_layout(fig)
        counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
        total = float(counts.sum()) or 1.0
        fracs = counts.values / total

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
            outside_legend(fig, ax, wedges, cats, legend_kws=legend_kws, title=groupby)
        ax.set(aspect="equal")
        if title:
            ax.set_title(title)
        save_or_show(fig, save, show=show)
        return ax

    splits = categorical_order(obs[split_by], split_order)
    n = len(splits)
    if ncols is None:
        ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    if figsize is None:
        figsize = (3.2 * ncols, 3.2 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    apply_layout(fig)
    axes = np.atleast_1d(axes).ravel()
    for i, sp in enumerate(splits):
        sub = obs[obs[split_by].astype(str) == sp]
        counts = sub[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
        ax = axes[i]
        total = float(counts.sum()) or 1.0
        fracs = counts.values / total

        wedges, texts, _ap = ax.pie(
            counts.values,
            labels=cats if show_labels else None,
            colors=[colors[c] for c in cats],
            autopct=_autopct if autopct else None,
            textprops=dict(fontsize=7),
            wedgeprops=dict(edgecolor="white", linewidth=0.5),
        )
        if show_labels:
            for t, f in zip(texts, fracs):
                if f < 0.02:
                    t.set_text("")
            if adjust_text:
                maybe_adjust_texts([t for t in texts if t.get_text()], ax=ax)
        ax.set_title(f"{sp} (n={int(total):,})", fontsize=9)
        ax.set(aspect="equal")
    for j in range(n, len(axes)):
        axes[j].axis("off")
    if show_legend:
        from matplotlib.patches import Patch

        handles = [Patch(color=colors[c], label=c) for c in cats]
        outside_legend(fig, axes[0], handles, legend_kws=legend_kws, title=groupby)
    if title:
        fig.suptitle(title)
    save_or_show(fig, save, show=show)
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
    figsize=None,
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    alpha: float = 0.85,
    show_legend: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
):
    """Stacked area chart of ``groupby`` composition across ``split_by``.

    For each ordered ``split_by`` category (x-axis) the ``groupby``
    proportions are stacked into filled bands, so the plot reads as the
    evolution of the composition across the ``split_by`` axis. By default
    the bands are column-normalised to fractions summing to 1.

    Parameters
    ----------
    adata : AnnData, pandas.DataFrame, or str
        Cell metadata source. When an :class:`~anndata.AnnData` is passed the
        function operates on ``adata.obs`` and can pick up colours from
        ``adata.uns[f"{groupby}_colors"]``; a DataFrame or a path to one is
        also accepted.
    groupby : str
        Column whose categories form the stacked bands.
    split_by : str
        Column defining the ordered x-axis positions (one stacked column
        per category).
    order : sequence, optional
        Explicit ordering of the ``groupby`` bands. Defaults to the natural
        categorical order of the column.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` categories along the x-axis.
        Defaults to the natural categorical order.
    normalize : bool, default True
        If ``True``, each column is normalised so the bands sum to 1
        (fractions); if ``False``, raw counts are stacked.
    palette : dict, str, or None, default None
        Colour mapping for the ``groupby`` bands. A dict maps
        category -> colour; a str is treated as a path to a palette
        workbook/sheet; ``None`` (default) resolves colours from the AnnData
        ``uns`` colours or a generated default.
    figsize : tuple of float, optional
        Figure size in inches, used only when a new figure is created
        (``ax`` is ``None``). When ``None`` (default) a 6 x 4 inch panel is
        widened for the outside legend and made tall enough to hold all of
        its entries.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw into. If ``None`` (default) a new figure and
        axes are created.
    save : str, optional
        If given, the figure is written to this path.
    show : bool, default False
        Whether to display the figure interactively.
    title : str, optional
        Axes title.
    alpha : float, default 0.85
        Opacity of the filled area bands.
    show_legend : bool, default True
        Whether to draw a legend of the ``groupby`` categories.
    legend_kws : dict, optional
        Extra keyword arguments forwarded to
        :meth:`matplotlib.axes.Axes.legend` (merged with defaults).

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the stacked area plot.
    """
    obs, ad = resolve_adata_obs(adata)
    splits = categorical_order(obs[split_by], split_order)
    cats = categorical_order(obs[groupby], order)
    ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
    ct = ct.reindex(index=cats, columns=splits, fill_value=0).astype(float)
    if normalize:
        col_sums = ct.sum(axis=0).replace(0, np.nan)
        ct = ct.div(col_sums, axis=1).fillna(0)

    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
    if figsize is None:
        figsize = _legend_aware_figsize((6, 4), cats, show_legend)
    own_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    apply_layout(fig, own=own_fig)

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
    ax.set_xticklabels(
        splits, rotation=-45, rotation_mode="anchor", ha="left", va="center"
    )
    ax.set_xlim(0, max(len(splits) - 1, 1))
    ax.set_ylim(0, 1 if normalize else None)
    ax.set_ylabel("Fraction" if normalize else "Count")
    ax.set_xlabel(split_by)
    ax.spines[["top", "right"]].set_visible(False)
    if show_legend:
        outside_legend(fig, ax, legend_kws=legend_kws, title=groupby)
    if title:
        ax.set_title(title)
    save_or_show(fig, save, show=show)
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
    """Dot plot where dot size = count and colour = column-fraction.

    Rows are ``groupby`` categories and columns are ``split_by`` categories.
    Each dot's area encodes the absolute cell count while its colour encodes
    the fraction of that ``split_by`` column made up by the ``groupby``
    category, making it easy to read both abundance and relative composition
    at once.

    Parameters
    ----------
    adata : AnnData, pandas.DataFrame, or str
        Cell metadata source. When an :class:`~anndata.AnnData` is passed the
        function operates on ``adata.obs``; a DataFrame or a path to one is
        also accepted.
    groupby : str
        Column whose categories form the plot rows.
    split_by : str
        Column whose categories form the plot columns; colour is the
        fraction within each of these columns.
    order : sequence, optional
        Explicit ordering of the ``groupby`` rows. Defaults to the natural
        categorical order of the column.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` columns. Defaults to the
        natural categorical order.
    palette : dict, str, or None, default None
        Accepted for API consistency with the other composition plots. The
        dot colours are driven by ``cmap`` (fraction values), so this is not
        used to colour the dots.
    cmap : str, default "viridis"
        Matplotlib colormap name used to encode the within-column fraction
        (mapped from 0 to the largest fraction in the table). A size legend
        with reference cell counts is drawn under the colour bar.
    size_max : float, default 250.0
        Marker area (in points squared) assigned to the largest count; all
        other dots are scaled linearly relative to it.
    figsize : tuple of float, optional
        Figure size in inches. When ``None`` (default) a size is derived
        from the number of rows and columns.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw into. If ``None`` (default) a new figure and
        axes are created.
    save : str, optional
        If given, the figure is written to this path.
    show : bool, default False
        Whether to display the figure interactively.
    title : str, optional
        Axes title.
    cbar_label : str, optional
        Label for the colour bar. Defaults to ``f"Fraction within {split_by}"``.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the dot plot.
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
        figsize = (
            max(4, 0.5 * len(splits) + 2) + 0.6 * label_extent(cats, 8),
            max(3, 0.35 * len(cats) + 1.5),
        )
    own_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    apply_layout(fig, own=own_fig)

    yy, xx = np.meshgrid(np.arange(len(cats)), np.arange(len(splits)), indexing="ij")
    max_count = max(counts.max(), 1)
    sizes = (counts / max_count) * size_max
    # Colour range follows the data: with many groups every fraction is
    # small and a fixed 0-1 range renders every dot in the same dark colour.
    vmax = float(frac.max()) if frac.size and frac.max() > 0 else 1.0
    sc = ax.scatter(
        xx.ravel(),
        yy.ravel(),
        s=sizes.ravel(),
        c=frac.ravel(),
        cmap=cmap,
        edgecolor="black",
        linewidth=0.3,
        vmin=0,
        vmax=vmax,
    )
    ax.set_xticks(np.arange(len(splits)))
    ax.set_xticklabels(
        splits, rotation=-45, rotation_mode="anchor", ha="left", va="center"
    )
    ax.set_yticks(np.arange(len(cats)))
    ax.set_yticklabels(cats)
    ax.set_xlim(-0.5, len(splits) - 0.5)
    ax.set_ylim(-0.5, len(cats) - 0.5)
    ax.invert_yaxis()
    ax.set_xlabel(split_by)
    ax.set_ylabel(groupby)
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.04, shrink=0.5, anchor=(0, 1))
    cbar.set_label(cbar_label or f"Fraction within {split_by}", fontsize=8)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=2, width=0.5, labelsize=7)
    # Size legend (cell counts) under the colorbar, outside the axes.
    handles, labels = sc.legend_elements(
        prop="sizes", num=4, func=lambda s: s / size_max * max_count,
        fmt="{x:,.0f}", color="grey", markeredgecolor="black", markeredgewidth=0.3,
    )
    if handles:
        ax.legend(
            handles, labels, title="Cells", loc="upper left",
            bbox_to_anchor=(1.02, 0.45), frameon=False, fontsize=7,
            title_fontsize=8, labelspacing=1.2, handletextpad=0.8, borderpad=0.6,
        )
    if title:
        ax.set_title(title)
    save_or_show(fig, save, show=show)
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
    figsize=None,
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    marker: str = "o",
    show_legend: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
):
    """Line trend of each ``groupby`` category across ordered ``split_by``.

    Plots one line per ``groupby`` category tracing its value across the
    ordered ``split_by`` axis. By default the values are column-normalised
    fractions, so the plot shows how the relative abundance of each category
    changes across conditions; set ``normalize=False`` for raw counts.

    Parameters
    ----------
    adata : AnnData, pandas.DataFrame, or str
        Cell metadata source. When an :class:`~anndata.AnnData` is passed the
        function operates on ``adata.obs`` and can pick up colours from
        ``adata.uns[f"{groupby}_colors"]``; a DataFrame or a path to one is
        also accepted.
    groupby : str
        Column whose categories are drawn as separate lines.
    split_by : str
        Column defining the ordered x-axis positions.
    order : sequence, optional
        Explicit ordering of the ``groupby`` lines. Defaults to the natural
        categorical order of the column.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` categories along the x-axis.
        Defaults to the natural categorical order.
    normalize : bool, default True
        If ``True``, each x position is column-normalised so the lines
        represent fractions; if ``False``, raw counts are plotted.
    palette : dict, str, or None, default None
        Colour mapping for the lines. A dict maps category -> colour; a str
        is treated as a path to a palette workbook/sheet; ``None`` (default)
        resolves colours from the AnnData ``uns`` colours or a generated
        default.
    figsize : tuple of float, optional
        Figure size in inches, used only when a new figure is created
        (``ax`` is ``None``). When ``None`` (default) a 6 x 4 inch panel is
        widened for the outside legend and made tall enough for its entries.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw into. If ``None`` (default) a new figure and
        axes are created.
    save : str, optional
        If given, the figure is written to this path.
    show : bool, default False
        Whether to display the figure interactively.
    title : str, optional
        Axes title.
    marker : str, default "o"
        Matplotlib marker style drawn at each data point along the lines.
    show_legend : bool, default True
        Whether to draw a legend of the ``groupby`` categories.
    legend_kws : dict, optional
        Extra keyword arguments forwarded to
        :meth:`matplotlib.axes.Axes.legend` (merged with defaults).

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the trend plot.
    """
    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    splits = categorical_order(obs[split_by], split_order)
    ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
    ct = ct.reindex(index=cats, columns=splits, fill_value=0).astype(float)
    if normalize:
        col_sums = ct.sum(axis=0).replace(0, np.nan)
        ct = ct.div(col_sums, axis=1).fillna(0)

    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
    if figsize is None:
        figsize = _legend_aware_figsize((6, 4), cats, show_legend)
    own_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    apply_layout(fig, own=own_fig)

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
    ax.set_xticklabels(
        splits, rotation=-45, rotation_mode="anchor", ha="left", va="center"
    )
    ax.set_xlim(-0.3, len(splits) - 0.7)
    ax.set_xlabel(split_by)
    ax.set_ylabel("Fraction" if normalize else "Count")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linewidth=0.4, alpha=0.5)
    if show_legend:
        outside_legend(fig, ax, legend_kws=legend_kws, title=groupby)
    if title:
        ax.set_title(title)
    save_or_show(fig, save, show=show)
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

    Parameters
    ----------
    adata : AnnData, pandas.DataFrame, or str
        Cell metadata source. When an :class:`~anndata.AnnData` is passed the
        function operates on ``adata.obs`` and can pick up colours from
        ``adata.uns[f"{groupby}_colors"]``; a DataFrame or a path to one is
        also accepted.
    groupby : str
        Column whose categories are stacked within each bar.
    split_by : str
        Column defining the individual bars (one bar per category).
    order : sequence, optional
        Explicit ordering of the ``groupby`` stack segments. Defaults to the
        natural categorical order of the column.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` bars along the x-axis.
        Defaults to the natural categorical order (may be overridden by
        ``sort_by``).
    normalize : bool, default True
        If ``True``, each bar is row-normalised so its segments sum to 1
        (proportions); if ``False``, raw counts are stacked.
    palette : dict, str, or None, default None
        Colour mapping for the ``groupby`` segments. A dict maps
        category -> colour; a str is treated as a path to a palette
        workbook/sheet; ``None`` (default) resolves colours from the AnnData
        ``uns`` colours or a generated default.
    figsize : tuple of float, optional
        Figure size in inches. When ``None`` (default) a size is derived
        from the number of bars and categories, plus room for the outside
        legend and the rotated bar labels.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw into. If ``None`` (default) a new figure and
        axes are created.
    save : str, optional
        If given, the figure is written to this path.
    show : bool, default False
        Whether to display the figure interactively (forwarded to the
        save/show helper).
    title : str, optional
        Axes title.
    rotation : float, default -45
        Rotation angle (degrees) of the x tick labels; the sign also
        controls their horizontal alignment.
    gap : float, default 0.05
        Fractional gap between adjacent bars (bar width is ``1 - gap``).
    edgecolor : str, default "black"
        Colour of the bar segment borders.
    linewidth : float, default 0.1
        Width of the bar segment borders.
    show_legend : bool, default True
        Whether to draw a legend of the ``groupby`` categories (placed
        outside the axes on the right).
    legend_kws : dict, optional
        Extra keyword arguments that override the default legend styling
        passed to :meth:`matplotlib.axes.Axes.legend`.
    sort_by : str, optional
        Name of a ``groupby`` category; when given, bars are sorted in
        ascending order of that category's value, overriding
        ``split_order``.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the stacked bar plot.
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
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
    if figsize is None:
        # Bars get a fixed width each; the legend column (longest label) and
        # the rotated x labels are added on top so they can't squeeze the
        # bars into slivers.
        width = max(2.5, len(splits) * 0.45) + 0.4 * label_extent(splits, 8, rotation)
        if show_legend:
            width += label_extent(cats, 8) + 0.4
        figsize = (width, max(3.5, min(len(cats) * 0.25 + 1.0, 8)))
    own_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    apply_layout(fig, own=own_fig)
    n = len(splits)
    bottom = np.zeros(n, dtype=float)
    x = np.arange(n, dtype=float)  # left edges
    bar_width = 1.0 - gap
    for c in cats:
        vals = ct[c].values
        ax.bar(
            x,
            vals,
            width=bar_width,
            bottom=bottom,
            align="edge",
            color=colors[c],
            edgecolor=edgecolor,
            linewidth=linewidth,
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
        ha="left" if rotation < 0 else ("right" if rotation else "center"),
        va="center" if rotation else "top",
    )
    ax.tick_params(
        axis="y",
        which="both",
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
    )
    ax.tick_params(axis="x", length=3, pad=2, top=False, labeltop=False)
    ax.xaxis.label.set_visible(False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(False)
    if show_legend:
        lkws = dict(
            loc="upper left",
            bbox_to_anchor=(1.0, 1.0),
            frameon=False,
            fontsize=8,
            ncol=1,
            borderpad=0.3,
            handlelength=1.0,
            handleheight=1.0,
            handletextpad=0.4,
            labelspacing=0.25,
            columnspacing=0.8,
            title=groupby,
            title_fontsize=9,
        )
        if legend_kws:
            lkws.update(legend_kws)
        ax.legend(**lkws)
    if title:
        ax.set_title(title)
    save_or_show(fig, save, show=show)
    return ax
