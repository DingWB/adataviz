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

    Supports 2 or 3 groups; for higher cardinality use :func:`upset_plot`.
    Requires the optional dependency :mod:`matplotlib_venn`.
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
    save_or_show(fig, save)
    return ax


def upset_plot(
    adata,
    groupby: str,
    set_by: str,
    order: Optional[Sequence] = None,
    min_subset_size: Optional[int] = None,
    figsize=(8, 4),
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    sort_by: str = "cardinality",
):
    """UpSet plot of ``set_by`` overlaps across ``groupby`` groups.

    Requires the optional dependency :mod:`upsetplot`.
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
    fig = plt.figure(figsize=figsize)
    kw = dict(sort_by=sort_by, show_counts=True)
    if min_subset_size is not None:
        kw["min_subset_size"] = min_subset_size
    upsetplot.UpSet(series, **kw).plot(fig=fig)
    if title:
        fig.suptitle(title)
    save_or_show(fig, save)
    return fig
