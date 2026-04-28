"""Pseudotime plot (violin)."""

from __future__ import annotations

import os
from typing import Any, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._utils import categorical_order, merge_legend_kws, save_or_show

__all__ = ["plot_pseudotime"]


def plot_pseudotime(
    pseudotime: Union[str, pd.DataFrame],
    groupby: str = "Age",
    y: str = "dpt_pseudotime",
    hue: Optional[str] = None,
    figsize=(5, 3.5),
    palette: Union[None, Mapping[str, str], str] = None,
    rotation: Optional[float] = None,
    ylabel: str = "Pseudotime",
    save: Optional[str] = None,
    show: bool = False,
    title: Optional[str] = None,
    legend_kws: Optional[Mapping[str, Any]] = None,
    sheet_name: Optional[str] = None,
):
    """Violin plot of pseudotime values grouped by a categorical variable.

    Parameters
    ----------
    pseudotime : path or DataFrame
        TSV file or DataFrame with at least the columns *y* and *groupby*.
    groupby : str
        Column to group cells along the X axis.
    y : str
        Pseudotime column.
    hue : str, optional
        Secondary grouping for split violins.
    palette : dict, path, or None
        Either a ``{group: hex}`` dict or an Excel palette path.
    legend_kws : dict, optional
        Forwarded to the legend (only used when ``hue`` is given).
    """
    import seaborn as sns

    if isinstance(pseudotime, str):
        data = pd.read_csv(os.path.expanduser(pseudotime), sep="\t", index_col=0)
    else:
        data = pseudotime.copy()
    if y in data.columns:
        data[y] = data[y].replace({np.inf: 1})
    cats = categorical_order(data[groupby])
    order = (
        data.groupby(groupby, observed=True)[y]
        .mean()
        .reindex(cats)
        .sort_values(ascending=True)
        .index.tolist()
    )
    hue_order = None
    if hue is not None:
        hue_order = (
            data.groupby(hue, observed=True)[y]
            .mean()
            .sort_values(ascending=True)
            .index.tolist()
        )

    if isinstance(palette, dict):
        color_palette = palette.copy()
    elif isinstance(palette, str) and os.path.exists(os.path.expanduser(palette)):
        sheet = sheet_name or groupby
        try:
            color_palette = pd.read_excel(
                os.path.expanduser(palette), sheet_name=sheet, index_col=0
            ).Hex.to_dict()
        except Exception:
            color_palette = None
    else:
        color_palette = None

    fig, ax = plt.subplots(figsize=figsize)
    sns.violinplot(
        data=data, x=groupby, y=y,
        order=order, hue=hue, hue_order=hue_order,
        palette=color_palette, density_norm="count", bw_adjust=0.5,
        inner=None, saturation=0.8, ax=ax, cut=0, linewidth=0.6,
    )
    ax.set_ylabel(ylabel)
    if rotation is not None:
        plt.setp(
            ax.get_xticklabels(),
            rotation=rotation,
            rotation_mode="anchor",
            ha="left" if rotation < 0 else "right",
        )
    if hue is not None:
        ax.legend(**merge_legend_kws(legend_kws, title=hue))
    if title:
        ax.set_title(title)
    save_or_show(fig, save, show=show)
    return ax
