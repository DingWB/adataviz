"""Flow-style plots: sankey (two-column flow) and chord (circular co-occurrence).

The sankey implementation is pure matplotlib (no extra dependencies). The
chord plot prefers :mod:`pycirclize` for nicer rendering and falls back
to a matplotlib polygon implementation when it is not installed.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._utils import (
    categorical_order,
    merge_legend_kws,
    resolve_adata_obs,
    resolve_palette,
    save_or_show,
)

__all__ = ["sankey_plot", "chord_plot"]


# ---------------------------------------------------------------------------
# Sankey
# ---------------------------------------------------------------------------


def _sankey_node_geometry(totals: np.ndarray, gap: float):
    """Return per-node ``(top, bottom, height)`` for a Sankey side.

    Nodes always span the full canvas: their heights are scaled to fit
    ``1 - gap * (n - 1)`` so that gaps between them never push the
    bottom node below 0 (this was the bug behind the empty-right-side
    artefact reported by users).
    """
    n = len(totals)
    if n == 0:
        return []
    available = max(0.01, 1.0 - gap * max(n - 1, 0))
    grand = max(totals.sum(), 1.0)
    sizes = totals / grand * available
    positions = []
    cur = 1.0
    for s in sizes:
        top = cur
        bottom = cur - s
        positions.append((top, bottom, s))
        cur = bottom - gap
    return positions


def sankey_plot(
    adata,
    left: str,
    right: str,
    left_order: Optional[Sequence] = None,
    right_order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=(8, 5),
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    gap: float = 0.015,
    flow_alpha: float = 0.55,
    label_fontsize: float = 9,
    show_legend: bool = False,
    legend_kws: Optional[Mapping[str, Any]] = None,
    min_flow_frac: float = 0.0,
    min_node_frac: float = 0.0,
):
    """Two-column Sankey/alluvial flow between ``left`` and ``right`` columns.

    Draws a Sankey (alluvial) diagram where each category on the ``left``
    column becomes a node on the left axis and each category on the
    ``right`` column becomes a node on the right axis. The width of every
    ribbon connecting a left node to a right node is proportional to the
    number of observations (rows of ``adata.obs``) that fall into that
    pair of categories, so the plot visualises how the ``left`` grouping
    maps onto the ``right`` grouping. Implemented as a pure matplotlib
    drawing: rectangles for nodes and smooth filled cubic-Bezier ribbons
    for flows.

    Parameters
    ----------
    adata : AnnData, DataFrame, or path
        Cell metadata source. Operates on ``adata.obs`` when given an
        AnnData; a DataFrame or path is resolved to a metadata table.
    left : str
        Column in ``adata.obs`` whose categories form the left-hand
        nodes of the diagram.
    right : str
        Column in ``adata.obs`` whose categories form the right-hand
        nodes of the diagram.
    left_order : sequence, optional
        Explicit ordering of the ``left`` categories from top to bottom.
        When ``None`` the natural categorical order is used.
    right_order : sequence, optional
        Explicit ordering of the ``right`` categories from top to bottom.
        When ``None`` the natural categorical order is used.
    palette : dict, str, or None, default None
        Colour specification for the ``left`` nodes and their ribbons.
        May be a ``{category: colour}`` mapping, a name/path resolved by
        the palette machinery, or ``None`` to fall back to
        ``adata.uns[<left>_colors]`` or a generated default palette. The
        ``right`` nodes are always coloured from the ``right`` column's
        own resolved palette.
    figsize : tuple of float, default (8, 5)
        Figure size in inches, used only when ``ax`` is None.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw into. When ``None`` a new figure and axes
        are created using ``figsize``.
    save : str, optional
        If given, path to write the figure to (delegated to
        ``save_or_show``); otherwise the figure is not saved.
    show : bool, default False
        Whether to display the figure interactively. Retained for API
        consistency; display is handled by ``save_or_show``.
    title : str, optional
        Title placed above the axes. No title is drawn when ``None``.
    gap : float, default 0.015
        Vertical gap (in axes-fraction units) inserted between adjacent
        nodes on each side. Node heights are scaled so the stack always
        fits within the canvas regardless of this gap.
    flow_alpha : float, default 0.55
        Opacity of the ribbons connecting the nodes.
    label_fontsize : float, default 9
        Base font size for node labels. Automatically shrunk when there
        are many nodes so labels do not overlap.
    show_legend : bool, default False
        When True, draw a legend of the ``left`` categories on the axes.
    legend_kws : dict, optional
        Forwarded to :meth:`matplotlib.axes.Axes.legend` when
        ``show_legend=True``.
    min_flow_frac : float, default 0.0
        Drop ribbons whose share of the grand total is below this
        fraction (e.g. ``0.005`` removes flows under 0.5%). Useful when
        the right column dwarfs a few large categories on the left.
    min_node_frac : float, default 0.0
        Drop nodes whose total share on either side is below this
        fraction. Smaller categories are pruned before layout so the
        remaining ones use the full canvas height.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the Sankey diagram.
    """
    obs, ad = resolve_adata_obs(adata)
    left_cats = categorical_order(obs[left], left_order)
    right_cats = categorical_order(obs[right], right_order)
    ct = pd.crosstab(obs[left].astype(str), obs[right].astype(str))
    ct = ct.reindex(index=left_cats, columns=right_cats, fill_value=0).astype(float)
    if ct.values.sum() == 0:
        raise ValueError("No data to plot in sankey_plot.")

    grand = float(ct.values.sum())
    # Optional filtering: drop tiny ribbons and tiny nodes for clarity.
    if min_flow_frac > 0:
        ct = ct.where(ct >= min_flow_frac * grand, 0.0)
    left_totals = ct.sum(axis=1).values
    right_totals = ct.sum(axis=0).values
    if min_node_frac > 0:
        keep_l = left_totals >= min_node_frac * grand
        keep_r = right_totals >= min_node_frac * grand
        left_cats = [c for c, k in zip(left_cats, keep_l) if k]
        right_cats = [c for c, k in zip(right_cats, keep_r) if k]
        ct = ct.loc[left_cats, right_cats]
        left_totals = ct.sum(axis=1).values
        right_totals = ct.sum(axis=0).values

    # Drop completely empty nodes — they would otherwise show as zero-height
    # rectangles with labels stacked on top of each other.
    keep_left = left_totals > 0
    keep_right = right_totals > 0
    if not keep_left.all():
        left_cats = [c for c, k in zip(left_cats, keep_left) if k]
        left_totals = left_totals[keep_left]
        ct = ct.loc[left_cats]
    if not keep_right.all():
        right_cats = [c for c, k in zip(right_cats, keep_right) if k]
        right_totals = right_totals[keep_right]
        ct = ct.loc[:, right_cats]

    left_pos = _sankey_node_geometry(left_totals, gap)
    right_pos = _sankey_node_geometry(right_totals, gap)

    colors = resolve_palette(
        palette, left_cats, sheet_name=left, adata=ad, groupby=left
    )
    right_colors = resolve_palette(
        None, right_cats, sheet_name=right, adata=ad, groupby=right
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Auto-shrink label fontsize when there are many nodes so labels don't overlap.
    fig_h_in = fig.get_size_inches()[1]
    max_n = max(len(left_cats), len(right_cats))
    auto_fs = max(4.0, min(label_fontsize, fig_h_in * 72 / max(max_n, 1) * 0.55))
    label_fontsize = float(auto_fs)

    # Draw nodes
    bar_w = 0.04
    for i, c in enumerate(left_cats):
        top, bottom, _ = left_pos[i]
        ax.add_patch(
            plt.Rectangle(
                (0, bottom), bar_w, top - bottom, facecolor=colors[c], edgecolor="none"
            )
        )
        ax.text(
            -0.01,
            (top + bottom) / 2,
            c,
            ha="right",
            va="center",
            fontsize=label_fontsize,
        )
    for j, c in enumerate(right_cats):
        top, bottom, _ = right_pos[j]
        ax.add_patch(
            plt.Rectangle(
                (1 - bar_w, bottom),
                bar_w,
                top - bottom,
                facecolor=right_colors[c],
                edgecolor="none",
            )
        )
        ax.text(
            1.01,
            (top + bottom) / 2,
            c,
            ha="left",
            va="center",
            fontsize=label_fontsize,
        )

    # Ribbons. Each ribbon's left thickness = v / left_total[i] * left_node_height
    # and right thickness = v / right_total[j] * right_node_height. This keeps
    # the ribbons aligned to their (shrunk) nodes even when ``gap`` reduces the
    # total height available for nodes.
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch

    left_cursor_top = [t for (t, _b, _h) in left_pos]
    right_cursor_top = [t for (t, _b, _h) in right_pos]
    left_node_h = np.array([h for (_t, _b, h) in left_pos])
    right_node_h = np.array([h for (_t, _b, h) in right_pos])

    for i, lc in enumerate(left_cats):
        lt = left_totals[i]
        if lt <= 0:
            continue
        for j, rc in enumerate(right_cats):
            v = ct.iat[i, j]
            if v <= 0:
                continue
            rt = right_totals[j]
            l_h = v / lt * left_node_h[i]
            r_h = v / rt * right_node_h[j]

            l_top = left_cursor_top[i]
            l_bot = l_top - l_h
            left_cursor_top[i] = l_bot
            r_top = right_cursor_top[j]
            r_bot = r_top - r_h
            right_cursor_top[j] = r_bot

            x0, x1 = bar_w, 1 - bar_w
            ctrl = (x0 + x1) / 2
            verts = [
                (x0, l_top),
                (ctrl, l_top),
                (ctrl, r_top),
                (x1, r_top),
                (x1, r_bot),
                (ctrl, r_bot),
                (ctrl, l_bot),
                (x0, l_bot),
                (x0, l_top),
            ]
            codes = [
                Path.MOVETO,
                Path.CURVE4,
                Path.CURVE4,
                Path.CURVE4,
                Path.LINETO,
                Path.CURVE4,
                Path.CURVE4,
                Path.CURVE4,
                Path.CLOSEPOLY,
            ]
            patch = PathPatch(
                Path(verts, codes),
                facecolor=colors[lc],
                edgecolor="none",
                alpha=flow_alpha,
            )
            ax.add_patch(patch)

    if show_legend:
        from matplotlib.patches import Patch

        handles = [Patch(color=colors[c], label=c) for c in left_cats]
        ax.legend(handles=handles, **merge_legend_kws(legend_kws, title=left))

    ax.set_xlim(-0.25, 1.25)
    ax.set_ylim(-0.05, 1.05)
    ax.axis("off")
    if title:
        ax.set_title(title)
    save_or_show(fig, save)
    return ax


# ---------------------------------------------------------------------------
# Chord
# ---------------------------------------------------------------------------


def chord_plot(
    adata,
    left: str,
    right: str,
    left_order: Optional[Sequence] = None,
    right_order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=(6, 6),
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    space: int = 2,
    adjust_text: bool = True,
    min_flow_frac: float = 0.0,
):
    """Circular chord diagram of ``left`` × ``right`` co-occurrence.

    Builds a symmetric co-occurrence matrix from the ``left`` and
    ``right`` categories and renders it as a chord (circular) diagram:
    every category becomes an arc (sector) around the ring, and the
    ribbons connecting two sectors are proportional to how often those
    two categories co-occur across observations in ``adata.obs``. This
    gives a compact circular view of the same relationships shown by
    :func:`sankey_plot`. Uses :mod:`pycirclize` when available (better
    labels and ribbons), otherwise falls back to a self-contained
    matplotlib renderer.

    Parameters
    ----------
    adata : AnnData, DataFrame, or path
        Cell metadata source. Operates on ``adata.obs`` when given an
        AnnData; a DataFrame or path is resolved to a metadata table.
    left : str
        Column in ``adata.obs`` providing the first set of categories.
    right : str
        Column in ``adata.obs`` providing the second set of categories.
    left_order : sequence, optional
        Explicit ordering of the ``left`` categories. When ``None`` the
        natural categorical order is used.
    right_order : sequence, optional
        Explicit ordering of the ``right`` categories. When ``None`` the
        natural categorical order is used.
    palette : dict, str, or None, default None
        Colour specification for the sectors. May be a
        ``{category: colour}`` mapping, a name/path resolved by the
        palette machinery, or ``None`` to fall back to
        ``adata.uns[<left>_colors]`` or a generated default palette.
    figsize : tuple of float, default (6, 6)
        Figure size in inches.
    save : str, optional
        If given, path to write the figure to (delegated to
        ``save_or_show``); otherwise the figure is not saved.
    show : bool, default False
        Whether to display the figure interactively. Retained for API
        consistency; display is handled by ``save_or_show``.
    title : str, optional
        Title placed above the figure. No title is drawn when ``None``.
    space : int, default 2
        Angular gap (in degrees) between adjacent sectors around the
        ring. Passed through to :class:`pycirclize.Circos` when
        available.
    adjust_text : bool, default True
        When True, nudge overlapping sector labels apart using
        ``adjustText`` (applies to both the pycirclize and matplotlib
        fallback renderers).
    min_flow_frac : float, default 0.0
        Drop ribbons whose share of the grand total is below this
        fraction; fully empty rows and columns are removed afterwards so
        zero-size sectors do not appear. Raises ``ValueError`` if the
        threshold removes every flow.

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the chord diagram.
    """
    obs, ad = resolve_adata_obs(adata)
    left_cats = categorical_order(obs[left], left_order)
    right_cats = categorical_order(obs[right], right_order)
    ct = pd.crosstab(obs[left].astype(str), obs[right].astype(str))
    ct = ct.reindex(index=left_cats, columns=right_cats, fill_value=0)
    if ct.values.sum() == 0:
        raise ValueError("No data to plot in chord_plot.")
    if min_flow_frac > 0:
        grand = float(ct.values.sum())
        ct = ct.where(ct >= min_flow_frac * grand, 0)
        if ct.values.sum() == 0:
            raise ValueError(
                f"min_flow_frac={min_flow_frac} removes all flows; lower the threshold."
            )
        # Drop fully-zero rows / columns so pycirclize does not divide by 0.
        keep_left = ct.sum(axis=1) > 0
        keep_right = ct.sum(axis=0) > 0
        ct = ct.loc[keep_left, keep_right]
        left_cats = [c for c in left_cats if c in ct.index]
        right_cats = [c for c in right_cats if c in ct.columns]

    all_labels = list(dict.fromkeys(list(left_cats) + list(right_cats)))
    matrix = pd.DataFrame(0.0, index=all_labels, columns=all_labels, dtype=float)
    for i, lc in enumerate(left_cats):
        for j, rc in enumerate(right_cats):
            v = float(ct.iat[i, j])
            if v == 0:
                continue
            matrix.at[lc, rc] = matrix.at[lc, rc] + v
            if lc != rc:
                matrix.at[rc, lc] = matrix.at[rc, lc] + v
    # Drop labels that ended up with zero total (avoids ZeroDivisionError
    # inside pycirclize.Circos when a sector has size 0).
    nz = matrix.sum(axis=1) > 0
    matrix = matrix.loc[nz, nz]
    all_labels = [l for l in all_labels if l in matrix.index]

    colors = resolve_palette(palette, all_labels, adata=ad, groupby=left)

    try:
        from pycirclize import Circos  # type: ignore
    except ImportError:
        return _chord_fallback(
            matrix,
            colors,
            figsize=figsize,
            save=save,
            show=show,
            title=title,
            adjust_text=adjust_text,
        )

    circos = Circos.initialize_from_matrix(
        matrix,
        space=space,
        cmap={k: colors[k] for k in all_labels},
        ticks_interval=None,
        label_kws=dict(size=8, r=105),
    )
    fig = circos.plotfig(figsize=figsize)
    if adjust_text:
        # Collect the sector-name Text artists pycirclize placed around the
        # ring and nudge them apart with adjustText.
        from ._utils import maybe_adjust_texts

        wanted = set(all_labels)
        texts = [t for ax in fig.axes for t in ax.texts if t.get_text() in wanted]
        if texts:
            maybe_adjust_texts(texts, ax=fig.axes[0])
    if title:
        fig.suptitle(title, y=0.98)
    fig.tight_layout()
    save_or_show(fig, save)
    return fig


def _chord_fallback(matrix, colors, figsize, save, show, title, adjust_text=True):
    """Minimal matplotlib chord fallback when pycirclize is unavailable."""
    from ._utils import maybe_adjust_texts

    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(aspect="equal"))
    labels = list(matrix.index)
    n = len(labels)
    totals = matrix.sum(axis=1).values
    grand = totals.sum()
    if grand <= 0:
        raise ValueError("Empty matrix.")

    gap = np.deg2rad(2)
    spans = (totals / grand) * (2 * np.pi - n * gap)
    starts = np.zeros(n)
    cur = 0.0
    for i in range(n):
        starts[i] = cur
        cur += spans[i] + gap
    ends = starts + spans

    R = 1.0
    R_outer = 1.05
    label_artists = []
    for i, lab in enumerate(labels):
        theta = np.linspace(starts[i], ends[i], 50)
        xs = np.concatenate([R * np.cos(theta), R_outer * np.cos(theta[::-1])])
        ys = np.concatenate([R * np.sin(theta), R_outer * np.sin(theta[::-1])])
        ax.fill(xs, ys, color=colors[lab], edgecolor="white", linewidth=0.5)
        mid = (starts[i] + ends[i]) / 2
        t = ax.text(
            1.12 * np.cos(mid),
            1.12 * np.sin(mid),
            lab,
            ha="center",
            va="center",
            fontsize=8,
            rotation=np.rad2deg(mid) - 90 if np.cos(mid) > 0 else np.rad2deg(mid) + 90,
        )
        label_artists.append(t)

    cursor_out = starts.copy()
    cursor_in = ends.copy()
    for i in range(n):
        for j in range(n):
            v = matrix.iat[i, j]
            if v <= 0 or i == j:
                continue
            frac_i = v / grand * (2 * np.pi - n * gap)
            a0 = cursor_out[i]
            a1 = cursor_out[i] + frac_i
            cursor_out[i] = a1
            b1 = cursor_in[j]
            b0 = b1 - frac_i
            cursor_in[j] = b0
            theta_a = np.linspace(a0, a1, 25)
            theta_b = np.linspace(b1, b0, 25)
            xs = np.concatenate(
                [R * np.cos(theta_a), [0.0], R * np.cos(theta_b), [0.0]]
            )
            ys = np.concatenate(
                [R * np.sin(theta_a), [0.0], R * np.sin(theta_b), [0.0]]
            )
            ax.fill(xs, ys, color=colors[labels[i]], alpha=0.35, edgecolor="none")
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.axis("off")
    if adjust_text:
        maybe_adjust_texts(label_artists, ax=ax)
    if title:
        ax.set_title(title)
    save_or_show(fig, save)
    return fig
