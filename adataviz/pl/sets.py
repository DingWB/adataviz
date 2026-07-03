"""Set-overlap plots: Venn (2/3 sets) and UpSet (any number).

These show how the values of one column (``set_by``) overlap across
groups defined by another column (``groupby``) in cell metadata.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt

from ._utils import (
    resolve_adata_obs,
    resolve_palette,
    save_or_show,
)

__all__ = ["venn_plot", "upset_plot"]


def _build_set_dict(adata, groupby, set_by):
    obs, _ = resolve_adata_obs(adata)
    out: dict = {}
    for g, sub in obs.groupby(groupby, observed=True):
        out[str(g)] = set(sub[set_by].dropna().astype(str).unique())
    return out


def venn_plot(
    adata,
    groupby: str,
    set_by: str,
    order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=(5, 5),
    ax=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
):
    """Venn diagram of ``set_by`` values present in each ``groupby`` group.

    For each group defined by ``groupby``, the unique (non-null) values in
    the ``set_by`` column are collected into a set, and the pairwise /
    triple-wise overlaps between those sets are drawn as a proportional
    Venn diagram. Supports 2 or 3 groups; for higher cardinality use
    :func:`upset_plot`. Requires the optional dependency
    :mod:`matplotlib_venn`.

    Parameters
    ----------
    adata : AnnData or path
        Data source. The observation metadata (``adata.obs`` or an
        equivalent table) must contain the ``groupby`` and ``set_by``
        columns; resolved via :func:`resolve_adata_obs`.
    groupby : str
        Metadata column whose categories define the sets (one circle per
        category). Must yield exactly 2 or 3 distinct groups after any
        ``order`` filtering.
    set_by : str
        Metadata column whose unique values populate each group's set;
        values are cast to ``str`` and null entries are dropped before the
        overlaps are computed.
    order : sequence, optional
        Explicit ordering / subset of the ``groupby`` categories to plot.
        Only entries that actually occur in the data are kept. When
        ``None`` (default), all groups are used in their natural order.
    palette : dict, str, or None, default None
        Colour mapping for the set circles. A ``{group: hex}`` dict, an
        Excel palette path/sheet name, or ``None`` to resolve colours from
        ``adata.uns`` / a default palette via :func:`resolve_palette`.
    figsize : tuple of float, default (5, 5)
        Figure size in inches; only used when ``ax`` is ``None`` and a new
        figure is created.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw into. When ``None`` (default) a new figure
        and axes are created.
    save : str, optional
        Path to write the figure to. When ``None`` (default) the figure is
        not saved.
    show : bool, default False
        Whether to display the figure interactively (forwarded to
        :func:`save_or_show`).
    title : str, optional
        Axes title. When ``None`` (default) no title is set.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the Venn diagram.
    """
    try:
        from matplotlib_venn import venn2, venn3  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "venn_plot requires `matplotlib-venn`. Install with "
            "`pip install matplotlib-venn`."
        ) from exc
    sets = _build_set_dict(adata, groupby, set_by)
    if order is not None:
        labels = [str(o) for o in order if str(o) in sets]
    else:
        labels = list(sets.keys())
    if len(labels) not in (2, 3):
        raise ValueError(
            f"venn_plot supports 2 or 3 sets; got {len(labels)} (use upset_plot)."
        )
    _obs, ad = resolve_adata_obs(adata)
    colors = resolve_palette(palette, labels, sheet_name=groupby, adata=ad, groupby=groupby)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    if len(labels) == 2:
        venn2(
            [sets[l] for l in labels],
            set_labels=labels,
            set_colors=[colors[l] for l in labels],
            ax=ax,
        )
    else:
        venn3(
            [sets[l] for l in labels],
            set_labels=labels,
            set_colors=[colors[l] for l in labels],
            ax=ax,
        )
    if title:
        ax.set_title(title)
    save_or_show(fig, save, show=show)
    return ax


def upset_plot(
    adata,
    groupby: str,
    set_by: str,
    order: Optional[Sequence] = None,
    min_subset_size: Optional[int] = None,
    figsize=None,
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    sort_by: str = "cardinality",
):
    """UpSet plot of ``set_by`` overlaps across ``groupby`` groups.

    Builds one set of ``set_by`` values per ``groupby`` category and
    visualises the size of every intersection as an UpSet plot (a matrix of
    set memberships paired with intersection-size bars). Scales to any
    number of groups, unlike :func:`venn_plot`. Requires the optional
    dependency :mod:`upsetplot`.

    Parameters
    ----------
    adata : AnnData or path
        Data source whose observation metadata contains ``groupby`` and
        ``set_by``; resolved via :func:`resolve_adata_obs`.
    groupby : str
        Metadata column whose categories define the sets (one column in the
        UpSet membership matrix per category).
    set_by : str
        Metadata column whose unique values populate each group's set;
        values are cast to ``str`` and null entries are dropped before the
        intersections are computed.
    order : sequence, optional
        Explicit ordering / subset of the ``groupby`` categories to include.
        Only entries present in the data are kept. When ``None`` (default),
        all groups are used in their natural order.
    min_subset_size : int, optional
        Minimum size an intersection must have to be shown; smaller
        intersections are hidden. When ``None`` (default) all non-empty
        intersections are drawn.
    figsize : tuple of float, optional
        Figure size in inches. When ``None`` (default) the width scales with
        the number of sets (``max(8.0, 0.7 * n_sets + 4.0)``) and the height
        is fixed at 5.0 so the matrix, bars and count labels stay readable.
    save : str, optional
        Path to write the figure to. When ``None`` (default) the figure is
        not saved.
    show : bool, default False
        Whether to display the figure interactively (forwarded to
        :func:`save_or_show`).
    title : str, optional
        Figure super-title. When ``None`` (default) no title is set.
    sort_by : str, default "cardinality"
        How the intersections are ordered along the plot, forwarded to
        :class:`upsetplot.UpSet` (e.g. ``"cardinality"`` to sort by
        intersection size or ``"degree"`` to sort by the number of
        participating sets).

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the UpSet plot.
    """
    try:
        import upsetplot  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "upset_plot requires `upsetplot`. Install with `pip install upsetplot`."
        ) from exc
    sets = _build_set_dict(adata, groupby, set_by)
    if order is not None:
        labels = [str(o) for o in order if str(o) in sets]
    else:
        labels = list(sets.keys())
    series = upsetplot.from_contents({l: sets[l] for l in labels})
    # Scale width with the number of sets so the matrix/bars and count
    # labels don't crowd together; give the combined layout enough height.
    if figsize is None:
        figsize = (max(8.0, 0.7 * len(labels) + 4.0), 5.0)
    fig = plt.figure(figsize=figsize)
    kw = dict(sort_by=sort_by, show_counts=True)
    if min_subset_size is not None:
        kw["min_subset_size"] = min_subset_size
    upsetplot.UpSet(series, **kw).plot(fig=fig)
    if title:
        fig.suptitle(title)
    save_or_show(fig, save, show=show)
    return fig
