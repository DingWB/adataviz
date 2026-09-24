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
    figsize=None,
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

    Draws one violin per ``groupby`` category showing the distribution of
    the pseudotime column ``y``. Groups are ordered along the X axis by
    their mean pseudotime (ascending), infinite pseudotime values are
    clamped to 1, and an optional ``hue`` produces split/grouped violins.

    Parameters
    ----------
    pseudotime : str or pandas.DataFrame
        Either a path to a tab-separated file (read with the first column as
        the index) or an in-memory DataFrame. Must contain at least the
        ``y`` and ``groupby`` columns (and ``hue`` when supplied).
    groupby : str, default "Age"
        Column used to group cells along the X axis; one violin is drawn per
        category, ordered by ascending mean of ``y``.
    y : str, default "dpt_pseudotime"
        Numeric pseudotime column plotted on the Y axis. Infinite values in
        this column are replaced with 1 before plotting.
    hue : str, optional
        Secondary categorical column for split/grouped violins. Its levels
        are ordered by ascending mean ``y``. When ``None`` (default) no hue
        split is applied and no legend is drawn.
    figsize : tuple of float, optional
        Figure size in inches. When ``None`` (default) the width scales with
        the number of violins (``max(5.0, 0.45 * n_groups + 1.5)``) and the
        height is fixed at 3.5 so densely populated axes stay legible.
    palette : dict, str, or None, default None
        Colour mapping for the violins. A ``{group: hex}`` dict is used
        directly; a path to an existing Excel palette file is read (using
        ``sheet_name`` or ``groupby`` as the sheet, with a ``Hex`` column);
        anything else (including ``None``) falls back to seaborn defaults.
    rotation : float, optional
        Rotation angle in degrees for the X tick labels. When ``None``
        (default) labels are auto-rotated to 45 degrees if there are more
        than 6 groups or any group name is longer than 6 characters, else
        left horizontal. Negative angles anchor the labels to the left.
    ylabel : str, default "Pseudotime"
        Label for the Y axis.
    save : str, optional
        Path to write the figure to. When ``None`` (default) the figure is
        not saved.
    show : bool, default False
        Whether to display the figure interactively (forwarded to
        :func:`save_or_show`).
    title : str, optional
        Axes title. When ``None`` (default) no title is set.
    legend_kws : dict, optional
        Extra keyword arguments merged (via :func:`merge_legend_kws`) into
        the legend call. Only used when ``hue`` is given; the legend title
        defaults to the ``hue`` column name.
    sheet_name : str, optional
        Sheet name to read when ``palette`` is an Excel file path. When
        ``None`` (default) the ``groupby`` column name is used as the sheet.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the pseudotime violin plot.
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

    # Scale width with the number of violins so they don't get crushed.
    if figsize is None:
        figsize = (max(5.0, 0.45 * len(order) + 1.5), 3.5)

    fig, ax = plt.subplots(figsize=figsize)
    sns.violinplot(
        data=data, x=groupby, y=y,
        order=order, hue=hue, hue_order=hue_order,
        palette=color_palette, density_norm="count", bw_adjust=0.5,
        inner=None, saturation=0.8, ax=ax, cut=0, linewidth=0.6,
    )
    ax.set_ylabel(ylabel)
    # Auto-rotate x tick labels when there are many or long category names,
    # so they don't overlap. Callers can still force a specific rotation.
    if rotation is None:
        maxlen = max((len(str(o)) for o in order), default=0)
        rotation = 45 if (len(order) > 6 or maxlen > 6) else 0
    if rotation:
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
