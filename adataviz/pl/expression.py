"""Per-gene expression plots that don't fit on the embedding.

- :func:`boxplot` -- multi-gene, multi-group box / violin / strip plot
  driven from a long-form pivot of ``adata.X``.
- :func:`stacked_violinplot` -- a clean seaborn-based stacked violin
  (one row per gene, one column-block per group).
- :func:`gene_dotplot` -- the rich PyComplexHeatmap dot heatmap that used
  to live as ``plot_genes`` (supports ``parent_col``, ``query_str``, etc.).
- :func:`gene_dotplot2D` -- one gene across two grouping variables
  (e.g. Subclass x Age), computed from a pseudobulk AnnData.
- :func:`get_genes_mean_frac` -- helper used by both ``gene_dotplot`` and
  the interactive variant.
"""

from __future__ import annotations

import os
import re
from typing import Any, Mapping, Optional, Sequence, Union

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._utils import (
    apply_layout,
    categorical_order,
    categorical_series,
    gene_positions,
    import_pch,
    label_extent,
    outside_legend,
    quiet_pch,
    resolve_adata_obs,
    resolve_palette,
    save_or_show,
)
from ..utils import normalize_mc_by_cell

__all__ = [
    "boxplot",
    "stacked_violinplot",
    "gene_dotplot",
    "gene_dotplot2D",
    "get_genes_mean_frac",
    "ridgeline",
]

# See complex_heatmap._DOT_LEGEND_KWS: keep the large reference dots of the
# size legend from overlapping their labels.
_DOT_LEGEND_KWS = dict(handlelength=3.0, handletextpad=0.8, labelspacing=0.9)


def _resolve_cutoff(df: pd.DataFrame, expression_cutoff):
    """Turn ``"p5"`` / ``"median"`` / ``"mean"`` into a numeric cutoff.

    Works on the flattened matrix without materialising a stacked Series
    (``df.stack()`` copies the data and is slow for many cells).
    """
    if not isinstance(expression_cutoff, str):
        return expression_cutoff
    values = df.to_numpy(dtype=float).ravel()
    values = values[~np.isnan(values)]
    if expression_cutoff == "median":
        return float(np.median(values))
    if expression_cutoff == "mean":
        return float(values.mean())
    return float(np.percentile(values, float(expression_cutoff.replace("p", ""))))


def _mean_frac_long(df: pd.DataFrame, groups, groupby: str, frac_matrix=None,
                    expression_cutoff=None, hypo: bool = False) -> pd.DataFrame:
    """Per-group mean and expressing fraction as a tidy ``[groupby, Gene, Mean, frac]``.

    ``frac_matrix`` (a cells x genes DataFrame already holding fractions, e.g. a
    pseudobulk ``frac`` layer) takes precedence; otherwise the fraction is the
    share of cells above ``expression_cutoff`` (or below 1 when ``hypo`` is set,
    for methylation-like data). Everything is a vectorised ``groupby().mean()``
    -- no per-group python lambdas and no tuple-keyed dict lookups.
    """
    mean_df = df.groupby(groups, observed=True).mean()
    if frac_matrix is not None:
        frac_df = frac_matrix.groupby(groups, observed=True).mean()
    elif hypo:
        frac_df = (df < 1).groupby(groups, observed=True).mean()
    else:
        frac_df = (df > expression_cutoff).groupby(groups, observed=True).mean()
    frac_df = frac_df.reindex(index=mean_df.index, columns=mean_df.columns)
    genes = list(mean_df.columns)
    cats = list(mean_df.index)
    return pd.DataFrame(
        {
            groupby: np.repeat(np.asarray(cats, dtype=object), len(genes)),
            "Gene": np.tile(np.asarray(genes, dtype=object), len(cats)),
            "Mean": mean_df.to_numpy(dtype=float).ravel(),
            "frac": frac_df.to_numpy(dtype=float).ravel(),
        }
    )


# ---------------------------------------------------------------------------
# Multi-gene / multi-group boxplot (or violin / strip)
# ---------------------------------------------------------------------------


def _gene_long_df(adata, genes, groupby, layer=None, use_raw=False):
    """Return a tidy ``[groupby, gene, value]`` DataFrame.

    ``groupby`` keeps the categorical order of ``adata.obs[groupby]`` (as an
    ordered ``Categorical`` of strings) so downstream plots follow the
    user's category order rather than alphabetical order.
    """
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("boxplot/stacked_violinplot require an AnnData.")
    src = adata.raw if use_raw and adata.raw is not None else adata
    # Column positions must come from the matrix we are about to slice: with
    # use_raw=True, adata.raw.var_names and adata.var_names can differ in both
    # order and length, so indices built from one silently select the wrong
    # genes out of the other.
    if layer is not None and layer in adata.layers:
        mat = adata.layers[layer]
        var = adata.var_names
    else:
        mat = src.X
        var = src.var_names
    X = mat[:, gene_positions(var, genes)]
    if hasattr(X, "toarray"):
        X = X.toarray()
    df = pd.DataFrame(np.asarray(X), index=adata.obs_names, columns=list(genes))
    df[groupby] = categorical_series(adata.obs[groupby]).values
    long = df.melt(id_vars=[groupby], var_name="gene", value_name="value")
    return long


def boxplot(
    adata,
    genes: Union[str, Sequence[str]],
    groupby: str,
    kind: str = "violin",
    layer: Optional[str] = None,
    use_raw: bool = False,
    order: Optional[Sequence] = None,
    gene_order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    figsize=None,
    ax=None,
    title: Optional[str] = None,
    show_strip: bool = False,
    strip_size: float = 0.5,
    show_legend: bool = True,
    legend_kws: Optional[Mapping[str, Any]] = None,
    rotation: float = -45,
    save: Optional[str] = None,
    show: bool = False,
    sharey: bool = False,
):
    """Multi-gene, multi-group box / violin / strip plot.

    Draws per-gene distributions of expression across the categories of
    ``groupby``. Expression values are pulled into a tidy long-form table
    (one row per cell x gene) from ``adata.X`` (or ``adata.raw.X`` /
    ``adata.layers[layer]``) and rendered with seaborn. A single gene
    produces one panel; multiple genes stack vertically, one panel per
    gene, sharing the group (X) axis.

    Parameters
    ----------
    adata : AnnData
        Source data. Must be an :class:`anndata.AnnData`; expression is
        read from it (see ``layer`` / ``use_raw``).
    genes : str or sequence of str
        Gene(s) to plot. Must exist in the resolved var space
        (``adata.raw.var_names`` when ``use_raw=True``, else
        ``adata.var_names``). A single string draws one panel; a sequence
        draws one stacked panel per gene.
    groupby : str
        Column in ``adata.obs`` used as the categorical X axis; each
        category becomes one violin/box/strip.
    kind : {"violin", "box", "strip"}, default ``"violin"``
        Plot style. ``"violin"`` uses :func:`seaborn.violinplot` with
        ``cut=0`` (distributions don't extend past the observed data
        range) and width-normalised densities; ``"box"`` uses
        :func:`seaborn.boxplot` with outliers hidden; ``"strip"`` uses
        :func:`seaborn.stripplot`.
    layer : str, optional
        Name of an ``adata.layers`` entry to read expression from. When
        given and present, it takes precedence over ``adata.X``.
    use_raw : bool, default False
        If True and ``adata.raw`` exists, read expression and gene names
        from ``adata.raw`` instead of ``adata``. Ignored when ``layer``
        is provided (the layer is always taken from ``adata``).
    order : sequence, optional
        Explicit ordering of ``groupby`` categories along the X axis.
        Categories not listed are dropped. When None, the natural
        categorical order of the data is used.
    gene_order : sequence, optional
        Explicit ordering of gene panels (top to bottom). Defaults to the
        order given in ``genes``.
    palette : dict, str, or None, optional
        Colour mapping for the ``groupby`` categories. A dict maps
        category -> colour; a string names a palette / colour sheet
        resolved via :func:`resolve_palette`. When None, colours are
        resolved from ``adata`` (e.g. stored ``*_colors``) or a default
        palette.
    figsize : tuple of float, optional
        Figure size in inches. When None it is derived from the number of
        categories and genes (wider for more groups, taller for more
        genes).
    ax : matplotlib.axes.Axes, optional
        Existing Axes to draw into. Only allowed when plotting a single
        gene; passing an Axes with multiple genes raises ``ValueError``.
        When None a new figure and Axes are created.
    title : str, optional
        Figure-level title drawn via ``fig.suptitle``.
    show_strip : bool, default False
        When ``kind`` is ``"violin"`` or ``"box"``, overlay a low-alpha
        black strip plot of the individual points on top.
    strip_size : float, default 0.5
        Marker size for strip points (both ``kind="strip"`` and the
        ``show_strip`` overlay).
    show_legend : bool, default True
        Whether to draw a category legend (only rendered on the first
        panel and only when ``legend_kws`` is provided).
    legend_kws : dict, optional
        Keyword arguments forwarded to :meth:`matplotlib.axes.Axes.legend`
        (merged with a ``title=groupby`` default). If None, no legend is
        drawn even when ``show_legend`` is True.
    rotation : float, default -45
        Rotation (degrees) applied to the X tick labels on the bottom
        panel. The sign controls horizontal alignment (left for negative,
        right otherwise).
    save : str, optional
        Path to write the figure to. When given the figure is saved via
        :func:`save_or_show`.
    show : bool, default False
        Whether to display the figure (``plt.show``) via
        :func:`save_or_show`.
    sharey : bool, default False
        When plotting multiple genes, whether the stacked panels share a
        common Y axis scale.

    Returns
    -------
    matplotlib.axes.Axes or list of matplotlib.axes.Axes
        The single Axes when one gene is plotted, otherwise the list of
        per-gene Axes.
    """
    import seaborn as sns

    if isinstance(genes, str):
        genes = [genes]
    long = _gene_long_df(adata, genes, groupby, layer=layer, use_raw=use_raw)
    cats = categorical_order(long[groupby], order)
    long = long[long[groupby].isin(cats)].copy()
    long[groupby] = pd.Categorical(long[groupby], categories=cats, ordered=True)

    _obs, ad = resolve_adata_obs(adata)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )

    if gene_order is None:
        gene_order = list(genes)

    n_genes = len(gene_order)
    if figsize is None:
        # Width grows with the number of groups; add room for rotated labels.
        figsize = (
            max(5, len(cats) * 0.5) + 0.4 * label_extent(cats, 8, rotation),
            2.5 * n_genes if n_genes > 1 else 3.5,
        )
    own_fig = ax is None
    if ax is None:
        if n_genes == 1:
            fig, axes = plt.subplots(1, 1, figsize=figsize)
            axes = [axes]
        else:
            fig, axes = plt.subplots(
                n_genes, 1, figsize=figsize, sharex=True, sharey=sharey
            )
            axes = list(np.atleast_1d(axes).ravel())
    else:
        if n_genes != 1:
            raise ValueError("Pass a single Axes only with one gene.")
        axes = [ax]
        fig = ax.figure
    # Constrained layout keeps the stacked panels tight while still making
    # room for the outside legend and the rotated labels of the last panel.
    apply_layout(fig, own=own_fig, h_pad=0.02, hspace=0.03)

    # seaborn >= 0.13 deprecates ``palette`` without ``hue``; colour by the
    # x variable itself (no dodging, no automatic duplicate legend).
    common = dict(
        x=groupby, y="value", order=cats, hue=groupby, hue_order=cats,
        palette=colors, legend=False, dodge=False,
    )
    for i, g in enumerate(gene_order):
        sub = long[long["gene"] == g]
        a = axes[i]
        if kind == "violin":
            sns.violinplot(
                data=sub, ax=a, cut=0, density_norm="width", linewidth=0.6,
                inner="quart", **common,
            )
        elif kind == "box":
            sns.boxplot(
                data=sub, ax=a, showfliers=False, linewidth=0.6, saturation=0.85,
                **common,
            )
        elif kind == "strip":
            sns.stripplot(
                data=sub, ax=a, size=strip_size, edgecolor="white", linewidth=0,
                **common,
            )
        else:
            raise ValueError(f"kind must be 'violin', 'box', or 'strip'; got {kind!r}")
        if show_strip and kind != "strip":
            sns.stripplot(
                data=sub,
                x=groupby,
                y="value",
                order=cats,
                color="black",
                ax=a,
                size=strip_size,
                alpha=0.4,
                jitter=0.3,
            )
        a.set_ylabel(g)
        a.spines[["top", "right"]].set_visible(False)
        if i < n_genes - 1:
            a.set_xlabel("")
            a.tick_params(labelbottom=False)
        else:
            a.set_xlabel(groupby)
            plt.setp(
                a.get_xticklabels(),
                rotation=rotation,
                rotation_mode="anchor",
                ha="left" if rotation < 0 else ("right" if rotation else "center"),
                va="center" if rotation else "top",
            )
        if show_legend and i == 0 and legend_kws is not None:
            from matplotlib.patches import Patch

            handles = [Patch(color=colors[c], label=c) for c in cats]
            # Figure-level "outside" legend: a legend taller than one panel
            # would otherwise force that panel's row to grow around it.
            outside_legend(fig, a, handles, legend_kws=legend_kws, title=groupby)
    if title:
        fig.suptitle(title)
    save_or_show(fig, save, show=show)
    return axes if n_genes > 1 else axes[0]


# ---------------------------------------------------------------------------
# Stacked violin (rewrite)
# ---------------------------------------------------------------------------


def stacked_violinplot(
    adata,
    genes: Sequence[str],
    groupby: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    standardize: Optional[str] = "gene",
    figsize=None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
    swap_axes: bool = False,
    rotation: Optional[float] = None,
    inner: Optional[str] = "quart",
):
    """Compact stacked violin plot.

    A clean seaborn-based stacked violin: each row (or each column when
    ``swap_axes=True``) is one gene and each tick along the shared axis is
    a ``groupby`` category. Optional per-gene or per-group min-max
    rescaling makes relative patterns comparable across genes. Requires no
    scanpy dependency.

    Parameters
    ----------
    adata : AnnData
        Source data. Expression is read from it (see ``layer`` /
        ``use_raw``).
    genes : sequence of str
        Genes to plot, one violin row (or column) each. Must exist in the
        resolved var space.
    groupby : str
        Column in ``adata.obs`` defining the categorical axis of groups.
    layer : str, optional
        Name of an ``adata.layers`` entry to read expression from; takes
        precedence over ``adata.X`` when present.
    use_raw : bool, default False
        If True and ``adata.raw`` exists, read expression and gene names
        from ``adata.raw``. Ignored when ``layer`` is supplied.
    order : sequence, optional
        Explicit ordering of ``groupby`` categories. Categories not listed
        are dropped; None uses the data's natural order.
    palette : dict, str, or None, optional
        Colour mapping for the ``groupby`` categories, resolved via
        :func:`resolve_palette` (dict, palette/sheet name, or None for a
        default).
    standardize : {"gene", "group", None}, default ``"gene"``
        Per-gene min-max rescaling to ``[0, 1]``. ``"gene"`` rescales each
        gene across all cells; ``"group"`` rescales each ``(group, gene)``
        combination independently; None leaves raw values unscaled.
    figsize : tuple of float, optional
        Figure size in inches. When None it is derived from the number of
        categories and genes and the ``swap_axes`` orientation.
    title : str, optional
        Figure-level title drawn via ``fig.suptitle``.
    save : str, optional
        Path to write the figure to (via :func:`save_or_show`).
    show : bool, default False
        Whether to display the figure (via :func:`save_or_show`).
    swap_axes : bool, default False
        When True, lay genes out as columns with horizontal violins
        instead of the default rows with vertical violins.
    rotation : float, optional
        Rotation (degrees) of the labels along the outer axis: the group
        tick labels (default ``-45``) or, with ``swap_axes=True``, the gene
        names under each narrow column (default ``90``, i.e. vertical, so
        neighbouring names cannot collide).
    inner : str, optional, default ``"quart"``
        ``inner`` argument forwarded to :func:`seaborn.violinplot`
        controlling the inner representation (e.g. ``"quart"`` for
        quartile lines, or None for nothing).

    Returns
    -------
    list of matplotlib.axes.Axes
        The per-gene Axes making up the stacked plot.
    """
    import seaborn as sns

    long = _gene_long_df(adata, genes, groupby, layer=layer, use_raw=use_raw)
    cats = categorical_order(long[groupby], order)
    long = long[long[groupby].isin(cats)].copy()
    long[groupby] = pd.Categorical(long[groupby], categories=cats, ordered=True)

    if standardize in ("gene", "group"):
        keys = "gene" if standardize == "gene" else [groupby, "gene"]
        grp = long.groupby(keys, observed=True)["value"]
        mn, mx = grp.transform("min"), grp.transform("max")
        rng = mx - mn
        # Constant groups (zero range) map to 0 instead of NaN.
        long["value"] = np.where(rng > 0, (long["value"] - mn) / rng.where(rng > 0, 1), 0.0)
    elif standardize is not None:
        raise ValueError("standardize must be 'gene', 'group' or None.")

    _obs, ad = resolve_adata_obs(adata)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
    common = dict(
        order=cats, hue=groupby, hue_order=cats, palette=colors, legend=False,
        dodge=False, cut=0, density_norm="width", linewidth=0.5, inner=inner,
    )

    n_genes = len(genes)
    if rotation is None:
        rotation = 90 if swap_axes else -45
    if swap_axes:
        if figsize is None:
            figsize = (
                max(4, n_genes * 0.5) + 0.5 * label_extent(cats, 8),
                # Vertical projection of the rotated gene names.
                max(3, len(cats) * 0.3) + 0.6 * label_extent(genes, 8, 90 - rotation),
            )
        fig, axes = plt.subplots(1, n_genes, figsize=figsize, sharey=True)
        apply_layout(fig, w_pad=0.02, wspace=0.03)
        axes = list(np.atleast_1d(axes).ravel())
        for i, g in enumerate(genes):
            sub = long[long["gene"] == g]
            sns.violinplot(data=sub, y=groupby, x="value", ax=axes[i], orient="h", **common)
            # Gene names hang below their column; anchored like rotated tick
            # labels (vertical by default) so neighbours don't collide.
            axes[i].set_xlabel(
                g, rotation=rotation, rotation_mode="anchor",
                ha=("left" if rotation < 0 else "right") if rotation % 180 else "center",
                va="center" if rotation % 180 else "top",
            )
            axes[i].set_ylabel("")
            axes[i].set_xticks([])
            axes[i].spines[["bottom", "top", "right"]].set_visible(False)
            axes[i].tick_params(axis="y", length=0)
        axes[0].set_ylabel(groupby)
    else:
        if figsize is None:
            figsize = (
                max(4, len(cats) * 0.4) + 0.5 * label_extent(genes, 8),
                max(2.5, n_genes * 0.5) + 0.4 * label_extent(cats, 8, rotation),
            )
        fig, axes = plt.subplots(n_genes, 1, figsize=figsize, sharex=True)
        apply_layout(fig, h_pad=0.02, hspace=0.03)
        axes = list(np.atleast_1d(axes).ravel())
        for i, g in enumerate(genes):
            sub = long[long["gene"] == g]
            sns.violinplot(data=sub, x=groupby, y="value", ax=axes[i], **common)
            axes[i].set_ylabel(g, rotation=0, ha="right", va="center", labelpad=8)
            axes[i].set_yticks([])
            axes[i].spines[["left", "top", "right"]].set_visible(False)
            if i < n_genes - 1:
                axes[i].set_xlabel("")
                axes[i].tick_params(labelbottom=False, length=0)
            else:
                axes[i].set_xlabel(groupby)
                plt.setp(
                    axes[i].get_xticklabels(),
                    rotation=rotation,
                    rotation_mode="anchor",
                    ha="left" if rotation < 0 else ("right" if rotation else "center"),
                    va="center" if rotation else "top",
                )
    if title:
        fig.suptitle(title)
    save_or_show(fig, save, show=show)
    return axes


# ---------------------------------------------------------------------------
# get_genes_mean_frac (helper)
# ---------------------------------------------------------------------------


def _resolve_obs_arg(adata, obs):
    if obs is None:
        return adata.obs.copy()
    if isinstance(obs, str):
        sep = "\t" if obs.endswith((".tsv", ".txt")) else ","
        return pd.read_csv(os.path.expanduser(obs), sep=sep, index_col=0)
    return obs.copy()


def get_genes_mean_frac(
    adata,
    genes: Sequence[str],
    groupby: str = "Subclass",
    obs=None,
    modality: str = "RNA",
    layer: str = "mean",
    use_raw: bool = False,
    expression_cutoff: Union[str, float] = "p5",
    normalize_per_cell: bool = True,
    clip_norm_value: float = 10,
    hypo_score: bool = False,
    query: Optional[str] = None,
):
    """Per-group mean expression and expressing-cell fraction.

    Computes, for every gene and every ``groupby`` category, the mean
    expression and the fraction of cells considered "expressing". This is
    the data-preparation helper behind :func:`gene_dotplot` and the
    interactive dot plot: it returns a tidy table rather than a figure.
    It transparently handles two input shapes -- ordinary cell-level
    AnnData and pre-aggregated pseudobulk AnnData (as produced by
    ``adataviz.tools.pseudobulk_stats``, detected by a ``"mean"`` layer),
    in which case stored ``mean`` and ``frac`` layers are used directly.

    Parameters
    ----------
    adata : AnnData or str
        Source data or a path to an ``.h5ad`` file (opened backed and
        closed after use). May be cell-level or pseudobulk.
    genes : sequence of str
        Genes to summarise. Genes absent from ``adata.var_names`` are
        reported and skipped.
    groupby : str, default ``"Subclass"``
        Grouping used for aggregation. For cell-level data it is a column
        in the (resolved) obs; for pseudobulk data it names the group
        label carried by ``obs_names`` and becomes the output column.
    obs : DataFrame, str, or None, optional
        Alternative cell metadata. A path is read as CSV/TSV (delimiter
        inferred from extension) with the first column as index; a
        DataFrame is used directly; None uses ``adata.obs``. Rows are
        intersected with ``adata.obs_names``.
    modality : str, default ``"RNA"``
        Data modality. For ``"RNA"``/``"ATAC"`` the expressing fraction
        uses an ``expression_cutoff``; other modalities (e.g. methylation)
        instead count cells below 1 as the fraction and may be normalised
        per cell.
    layer : str, default ``"mean"``
        For pseudobulk input, the layer used as the mean-expression
        matrix (paired with the ``"frac"`` layer for fractions).
    use_raw : bool, default False
        For cell-level input, if True and ``adata.raw`` exists, use raw
        counts as the expression matrix.
    expression_cutoff : str or float, default ``"p5"``
        Threshold above which a cell counts as expressing (RNA/ATAC). A
        float is used directly; a string is resolved from the data as a
        percentile (``"p5"`` -> 5th percentile), ``"median"``, or
        ``"mean"`` of all non-aggregated values.
    normalize_per_cell : bool, default True
        For non-RNA/ATAC modalities, whether to normalise each cell before
        aggregating (via :func:`normalize_mc_by_cell`).
    clip_norm_value : float, default 10
        Clip value passed to :func:`normalize_mc_by_cell` when
        normalising non-RNA/ATAC data.
    hypo_score : bool, default False
        Passed to :func:`normalize_mc_by_cell`; controls hypo-methylation
        scoring for methylation-like modalities.
    query : str, optional
        Pandas-style query applied to ``adata.obs`` (or the resolved
        ``obs``) to subset cells/groups *before* loading expression into
        memory -- useful with backed AnnData to avoid materialising the
        full matrix. For pseudobulk input the group label is exposed as a
        queryable column. Raises if the query matches no rows.

    Returns
    -------
    pandas.DataFrame
        Long-form table with columns ``[groupby, "Gene", "Mean", "frac"]``,
        where ``"Mean"`` is the per-group mean expression and ``"frac"``
        the fraction of expressing cells for each gene/group.
    """
    if isinstance(adata, str):
        adata = anndata.read_h5ad(os.path.expanduser(adata), backed="r")
    keep = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in keep]
    if missing:
        print(f"genes not found in adata: {missing}")

    # Pseudobulk detection. ``pseudobulk_stats`` (in adataviz.tools) emits
    # an adata where:
    #   * obs.index = unique values of the aggregation groupby (e.g.
    #     subclass labels), obs.index.name = that groupby.
    #   * obs.columns = ['cell_count'] only (no Subclass / Group column).
    #   * layers carry ``min/q25/q50/q75/max/sum/std/frac/mean``.
    # Cell-level obs has no meaning here; the groups are already
    # aggregated.
    is_pseudobulk = "mean" in adata.layers

    if query is not None:
        if is_pseudobulk:
            # Make the group label queryable as a column under any of:
            # the user's ``groupby`` name and the original
            # ``obs.index.name``. This lets ``query="Subclass=='X'"``
            # work even though pseudobulk_stats does not store
            # ``Subclass`` as an obs column.
            aug = adata.obs.copy()
            label_cols = {groupby, aug.index.name} - {None}
            for col in label_cols:
                if col not in aug.columns:
                    aug[col] = aug.index.astype(str)
            kept_cells = aug.query(query).index
        else:
            # Cell-level adata.
            if obs is not None:
                obs = _resolve_obs_arg(adata, obs).query(query)
                kept_cells = obs.index.intersection(adata.obs_names)
            else:
                kept_cells = adata.obs.query(query).index
        if len(kept_cells) == 0:
            raise ValueError(
                f"query={query!r} matched 0 rows in adata.obs "
                f"(pseudobulk={is_pseudobulk}). "
                f"Available labels sample: {list(adata.obs_names[:5])}."
            )
        # h5py forbids fancy indexing on both axes simultaneously, so
        # subset rows on the backed object first, materialize, then
        # subset columns in memory.
        use = adata[kept_cells].to_memory()[:, keep].copy()
    else:
        use = adata[:, keep].to_memory()
    if hasattr(adata, "isbacked") and adata.isbacked:
        adata.file.close()

    if is_pseudobulk or "mean" in use.layers:
        # Pseudobulk path. ``use.obs_names`` IS the list of group labels
        # (no need to consult ``use.obs`` for a groupby column — it
        # doesn't carry one). Build the long table straight from the two
        # aligned matrices so genes whose mean is NaN for a queried single
        # group are kept (``stack()`` would drop them, and
        # ``stack(dropna=False)`` is rejected by pandas >= 3).
        mean_df = use.to_df(layer=layer)
        frac_df = use.to_df(layer="frac").reindex(
            index=mean_df.index, columns=mean_df.columns
        )
        labels = list(mean_df.index.astype(str))
        gene_names = list(mean_df.columns)
        plot_data = pd.DataFrame(
            {
                groupby: np.repeat(np.asarray(labels, dtype=object), len(gene_names)),
                "Gene": np.tile(np.asarray(gene_names, dtype=object), len(labels)),
                "Mean": mean_df.to_numpy(dtype=float).ravel(),
                "frac": frac_df.to_numpy(dtype=float).ravel(),
            }
        )
        if plot_data.empty:
            raise ValueError(
                "Pseudobulk plot_data is empty — check that `keep` "
                "(intersection of `genes` with adata.var_names) is "
                "non-empty and that the queried group(s) exist."
            )
        return plot_data

    if use_raw and use.raw is not None:
        raw = use.raw[:, use.var_names].to_adata()
        use.X = raw[use.obs_names, use.var_names].X.copy()
    obs = _resolve_obs_arg(use, obs)
    overlap = use.obs_names[use.obs_names.isin(obs.index)]
    obs = obs.loc[overlap]
    use = use[overlap, :]

    if modality not in ("RNA", "ATAC") and normalize_per_cell:
        use = normalize_mc_by_cell(
            use_adata=use,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
            hypo_score=hypo_score,
        )

    df = use.to_df()
    if modality in ("RNA", "ATAC"):
        expression_cutoff = _resolve_cutoff(df, expression_cutoff)
    groups = categorical_series(obs.loc[df.index, groupby])
    plot_data = _mean_frac_long(
        df,
        groups,
        groupby,
        frac_matrix=use.to_df(layer="frac") if "frac" in use.layers else None,
        expression_cutoff=expression_cutoff,
        hypo=modality not in ("RNA", "ATAC"),
    )
    return plot_data


# ---------------------------------------------------------------------------
# gene_dotplot (rich, with parent_col / hierarchical annotation)
# ---------------------------------------------------------------------------


def _strip_cbar_white_lines(cm):
    """Remove the spurious white outline at the bottom of PyComplexHeatmap colorbars.

    PyComplexHeatmap draws its colorbars inside a child Axes whose spines
    sometimes render as a thin white line at the bottom under certain
    matplotlib + colorbar combinations. Hide all spines on each cax.
    """
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


def gene_dotplot(
    adata,
    genes: Sequence[str],
    groupby: Union[str, Sequence[str]] = "Subclass",
    parent_col: Optional[str] = None,
    obs=None,
    query_str: Optional[str] = None,
    modality: str = "RNA",
    use_raw: bool = True,
    expression_cutoff: Union[str, float] = "p5",
    cell_type_order: Optional[Sequence] = None,
    gene_order: Optional[Sequence] = None,
    row_cluster: bool = False,
    col_cluster: bool = False,
    cmap: str = "Reds",
    palette=None,
    legend_kws: Optional[Mapping[str, Any]] = None,
    normalize_per_cell: bool = True,
    clip_norm_value: float = 10,
    hypo_score: bool = False,
    figsize=None,
    marker: str = "o",
    plot_kws: Optional[Mapping[str, Any]] = None,
    transpose: bool = False,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
    outline=True,outline_kws: Optional[Mapping[str, Any]] = None
):
    """Rich PyComplexHeatmap dot heatmap with hierarchical annotations.

    Draws a dot heatmap where dot colour encodes mean expression and dot
    size encodes the fraction of expressing cells, for each gene across
    the categories of ``groupby``. Compared to a plain dot plot, this
    function adds an optional two-level ``parent_col`` annotation,
    ``query_str`` pre-filtering, palette-sheet colour resolution, and full
    ``plot_kws`` passthrough to ``DotClustermapPlotter`` -- aimed at
    figure-quality marker panels.

    Parameters
    ----------
    adata : AnnData or str
        Source data or a path to an ``.h5ad`` file (opened backed and
        closed after use).
    genes : sequence of str
        Genes to display (one heatmap row, or column when ``transpose``).
        Genes absent from ``adata.var_names`` are reported and skipped.
    groupby : str or sequence of str, default ``"Subclass"``
        Grouping used for the group axis. A sequence is joined with
        ``"+"`` into a single composite grouping (and each part gets its
        own annotation track).
    parent_col : str, optional
        A coarser ``adata.obs`` column used as a second-level (parent)
        annotation track and to order the groups hierarchically. When None
        only the ``groupby`` annotation is drawn.
    obs : DataFrame, str, or None, optional
        Alternative cell metadata (DataFrame or path); None uses
        ``adata.obs``. Rows are intersected with ``adata.obs_names``.
    query_str : str, optional
        Pandas-style query applied to ``obs`` to subset cells before
        aggregation. Also used as the default ``title`` when set.
    modality : str, default ``"RNA"``
        Data modality; controls the expressing-fraction computation and
        whether per-cell normalisation is applied (see ``get_genes_mean_frac``).
    use_raw : bool, default True
        If True and ``adata.raw`` exists, use raw counts as the expression
        matrix.
    expression_cutoff : str or float, default ``"p5"``
        Threshold above which a cell counts as expressing (RNA/ATAC). A
        float is used directly; a string is resolved as a percentile
        (``"p5"``), ``"median"``, or ``"mean"``.
    cell_type_order : sequence, optional
        Explicit ordering of the group categories along the group axis;
        entries not present in the data are ignored.
    gene_order : sequence, optional
        Explicit ordering of genes along the gene axis. Defaults to the
        kept-gene order.
    row_cluster : bool, default False
        Whether to hierarchically cluster heatmap rows.
    col_cluster : bool, default False
        Whether to hierarchically cluster heatmap columns.
    cmap : str, default ``"Reds"``
        Matplotlib colormap name mapping mean expression to dot colour.
    palette : dict, str, or None, optional
        Colours for the annotation tracks. A dict maps category ->
        colour for ``groupby``; a path to an ``.xlsx`` colour sheet is
        read per sheet (``groupby``, its ``"+"`` parts, and ``parent_col``);
        None resolves default palettes via :func:`resolve_palette`.
    legend_kws : dict, optional
        Colorbar legend options forwarded as ``cmap_legend_kws`` to
        ``DotClustermapPlotter``. Defaults to
        ``dict(extendfrac=0.1, extend="both", label="Mean")``.
    normalize_per_cell : bool, default True
        For non-RNA/ATAC modalities, whether to normalise each cell before
        aggregating (via :func:`normalize_mc_by_cell`).
    clip_norm_value : float, default 10
        Clip value passed to :func:`normalize_mc_by_cell`.
    hypo_score : bool, default False
        Passed to :func:`normalize_mc_by_cell` for methylation-like data.
    figsize : tuple of float, optional
        Figure size in inches. When None it is auto-sized from the number
        of genes and groups so dots and tick labels don't overlap.
    marker : str, default ``"o"``
        Marker shape used for the dots (and the dot-size legend).
    plot_kws : dict, optional
        Extra keyword arguments forwarded to
        :class:`PyComplexHeatmap.DotClustermapPlotter`. Sensible defaults
        (clustering, tick-label sides/rotation, legend placement, spines,
        etc.) are filled in via ``setdefault`` and can be overridden here.
    transpose : bool, default False
        When True, place genes on the X axis and groups on the Y axis
        (annotation moves to the left) instead of the default layout.
    title : str, optional
        Figure title. Defaults to ``query_str`` (or ``groupby`` when no
        query is given).
    save : str, optional
        Path to write the figure to (via :func:`save_or_show`).
    show : bool, default False
        Whether to display the figure (via :func:`save_or_show`).

    Returns
    -------
    PyComplexHeatmap.DotClustermapPlotter
        The fitted clustermap plotter object (its axes are accessible via
        ``.heatmap_axes`` etc.).
    """
    pch = import_pch()
    HeatmapAnnotation, anno_simple, DotClustermapPlotter = (
        pch.HeatmapAnnotation,
        pch.anno_simple,
        pch.DotClustermapPlotter,
    )

    plot_kws = dict(plot_kws or {})
    if legend_kws is None:
        legend_kws = dict(extendfrac=0.1, extend="both", label="Mean")

    raw = (
        anndata.read_h5ad(os.path.expanduser(adata), backed="r")
        if isinstance(adata, str)
        else adata
    )
    keep = [g for g in genes if g in raw.var_names]
    missing = [g for g in genes if g not in keep]
    if missing:
        print(f"genes not found in adata: {missing}")
    # Reading from a backed h5ad requires positional indices in increasing
    # order (h5py limitation); select genes in var order, then restore the
    # requested order in memory.
    keep_set = set(keep)
    keep_sorted = [g for g in raw.var_names if g in keep_set]
    use = raw[:, keep_sorted].to_memory()
    use = use[:, keep].copy()
    if use_raw and use.raw is not None:
        rraw = use.raw[:, use.var_names].to_adata()
        use.X = rraw[use.obs_names, use.var_names].X.copy()
    if isinstance(adata, str) and raw.isbacked:
        raw.file.close()

    obs = _resolve_obs_arg(use, obs)
    if query_str is not None:
        obs = obs.query(query_str)
    overlap = use.obs_names[use.obs_names.isin(obs.index)]
    obs = obs.loc[overlap]
    # ``.copy()`` materialises the view so subsequent ``use.obs[...] = ...``
    # assignments don't emit ``ImplicitModificationWarning``.
    use = use[overlap, :].copy()

    if isinstance(groupby, list):
        joined = "+".join(groupby)
        obs[joined] = obs.loc[:, groupby].apply(
            lambda x: "+".join(x.astype(str)), axis=1
        )
        groupby = joined
    use.obs[groupby] = obs.loc[use.obs_names, groupby].tolist()
    if title is None:
        title = query_str or groupby

    if parent_col is not None and parent_col not in use.obs.columns:
        use.obs[parent_col] = obs.loc[use.obs_names, parent_col].tolist()
    if modality not in ("RNA", "ATAC") and normalize_per_cell:
        use = normalize_mc_by_cell(
            use_adata=use,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
            hypo_score=hypo_score,
        )

    color_palette: dict = {}
    if isinstance(palette, dict):
        color_palette[groupby] = palette.copy()
    elif isinstance(palette, str) and os.path.exists(os.path.expanduser(palette)):
        D = pd.read_excel(os.path.expanduser(palette), sheet_name=None, index_col=0)
        if groupby in D:
            color_palette[groupby] = D[groupby].Hex.to_dict()
        else:
            for sub in groupby.split("+"):
                if sub in D:
                    color_palette[sub] = D[sub].Hex.to_dict()
        if parent_col is not None and parent_col in D:
            color_palette[parent_col] = D[parent_col].Hex.to_dict()
    elif palette is None:
        cats = sorted(use.obs[groupby].astype(str).unique())
        color_palette[groupby] = resolve_palette(
            None,
            cats,
            sheet_name=groupby,
            adata=use if isinstance(use, anndata.AnnData) else None,
            groupby=groupby,
        )
        if parent_col is not None:
            pcats = sorted(use.obs[parent_col].astype(str).unique())
            color_palette[parent_col] = resolve_palette(
                None,
                pcats,
                adata=use,
                groupby=parent_col,
            )

    def _colors_for(key: str, cats) -> dict:
        """Colors for one annotation column, filling whatever the palette lacks.

        A "+"-joined ``groupby`` is annotated one sub-column at a time, so each
        sub-column needs its own palette even though the caller supplied one only
        for the composite key (or none at all).
        """
        cats = [str(c) for c in cats]
        have = color_palette.get(key) or {}
        if not all(c in have for c in cats):
            generated = resolve_palette(
                palette,
                cats,
                sheet_name=key,
                adata=use if isinstance(use, anndata.AnnData) else None,
                groupby=key,
            )
            have = {**generated, **have}
            color_palette[key] = have
        return {c: have[c] for c in cats}

    df = use.to_df()
    if modality in ("RNA", "ATAC"):
        expression_cutoff = _resolve_cutoff(df, expression_cutoff)
    # Groups follow the categorical order of the obs column when it has one
    # (falling back to sorted labels), both for the aggregation and for the
    # column order of the plot.
    obs_groups = obs.loc[use.obs_names, groupby]
    group_order = categorical_order(obs_groups)
    groups = categorical_series(obs_groups, group_order)
    groups.index = df.index
    plot_data = _mean_frac_long(
        df,
        groups,
        groupby,
        frac_matrix=use.to_df(layer="frac") if "frac" in use.layers else None,
        expression_cutoff=expression_cutoff,
        hypo=modality not in ("RNA", "ATAC"),
    )

    df_cols = pd.DataFrame({groupby: group_order}, index=group_order)
    if parent_col is not None:
        group2parent = (
            use.obs.loc[:, [groupby, parent_col]]
            .astype(str)
            .drop_duplicates(groupby)
            .set_index(groupby)[parent_col]
        )
        parent_src = obs if parent_col in obs.columns else use.obs
        parent_order = categorical_order(parent_src.loc[use.obs_names, parent_col])
        df_cols[parent_col] = df_cols[groupby].map(group2parent)
        df_cols = df_cols[df_cols[parent_col].isin(parent_order)]
        # Stable sort keeps the group order inside each parent block.
        rank = {p: i for i, p in enumerate(parent_order)}
        df_cols = df_cols.iloc[np.argsort(df_cols[parent_col].map(rank).values, kind="stable")]
    if cell_type_order is not None:
        df_cols = df_cols.loc[[c for c in cell_type_order if c in df_cols.index]]
    df_cols = df_cols.dropna()

    axis = 1 if not transpose else 0

    if "+" in groupby:
        anno_dict = {}
        for sub in groupby.split("+"):
            df_cols[sub] = df_cols[groupby].apply(
                lambda x, _s=sub: x.split("+")[groupby.split("+").index(_s)]
            )
            sub_colors = _colors_for(sub, df_cols[sub].unique())
            anno_dict[sub] = anno_simple(
                df_cols[sub],
                colors=sub_colors,
                add_text=False,
                legend=False,
                height=3,
                label=sub,
            )
        if parent_col is not None:
            pcolors = _colors_for(parent_col, df_cols[parent_col].unique())
            anno_dict[parent_col] = anno_simple(
                df_cols[parent_col],
                colors=pcolors,
                add_text=False,
                legend=True,
                height=3,
                label=parent_col,
            )
        col_ha = HeatmapAnnotation(**anno_dict, axis=axis, verbose=0)
    else:
        gcolors = _colors_for(groupby, df_cols[groupby].unique())
        kw = dict(
            axis=axis,
            group=anno_simple(
                df_cols[groupby],
                colors=gcolors,
                add_text=False,
                legend=False,
                height=3,
                label=groupby,
            ),
        )
        if parent_col is not None:
            pcolors = _colors_for(parent_col, df_cols[parent_col].unique())
            kw["parent"] = anno_simple(
                df_cols[parent_col],
                colors=pcolors,
                add_text=False,
                legend=True,
                height=3,
                label=parent_col,
            )
        col_ha = HeatmapAnnotation(**kw, verbose=0)

    # Always show tick labels; pick a side that does not collide with col_ha.
    # When col_ha is on top (axis=1, default) we keep group names at the
    # bottom of the heatmap and gene names on the right - both sides of the
    # heatmap stay clean and labels are clearly readable.
    if not transpose:
        top_anno, left_anno = col_ha, None
        x, y = groupby, "Gene"
        x_order = df_cols.index.tolist()
        y_order = gene_order if gene_order is not None else keep
        col_side, row_side = "bottom", "left"
    else:
        top_anno, left_anno = None, col_ha
        x, y = "Gene", groupby
        x_order = gene_order if gene_order is not None else keep
        y_order = df_cols.index.tolist()
        col_side, row_side = "bottom", "right"

    # Auto-size the figure from the grid dimensions so tick labels and dots
    # don't overlap when there are many genes/groups. (A fixed size crushes
    # labels into each other; bbox_inches="tight" can't fix that.)
    if figsize is None:
        n_groups = len(df_cols.index)
        n_genes = len(keep)
        n_x, n_y = (n_groups, n_genes) if not transpose else (n_genes, n_groups)
        figsize = (max(6.0, 0.30 * n_x + 3.0), max(3.0, 0.24 * n_y + 2.0))

    defaults = dict(
        marker=marker,
        grid=None,
        dot_legend_marker=marker,
        cmap_legend_kws=legend_kws,
        dot_legend_kws=dict(_DOT_LEGEND_KWS, title="Fraction"),
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        row_cluster_method="ward",
        row_cluster_metric="euclidean",
        col_cluster_method="ward",
        col_cluster_metric="euclidean",
        col_names_side=col_side,
        row_names_side=row_side,
        show_rownames=True,
        show_colnames=True,
        row_dendrogram=False,
        xticklabels_kws=dict(labelrotation=-45, labelsize=8, bottom=True),
        yticklabels_kws=dict(labelsize=8, right=True),
        legend_hpad=4,
        legend_side="right",
        legend_anchor="ax_heatmap",
        spines=False,
    )
    for k, v in defaults.items():
        plot_kws.setdefault(k, v)
    if title:
        plot_kws.setdefault("xlabel", title)
        plot_kws.setdefault("xlabel_side", 'top')
    fig = plt.figure(figsize=figsize)
    with quiet_pch(fig):
        cm = DotClustermapPlotter(
            data=plot_data,
            top_annotation=top_anno,
            left_annotation=left_anno,
            # Explicit orders only apply to un-clustered axes.
            x_order=None if plot_kws["col_cluster"] else x_order,
            y_order=None if plot_kws["row_cluster"] else y_order,
            x=x,
            y=y,
            value="Mean",
            c="Mean",
            s="frac",
            cmap=cmap,
            verbose=0,
            **plot_kws,
        )
    _style_dot_axes(cm, outline, outline_kws)
    with quiet_pch(fig):
        save_or_show(fig, save, show=show)
    return cm


def _style_dot_axes(cm, outline: bool, outline_kws: Optional[Mapping[str, Any]]):
    """Dashed minor grid behind the dots plus an optional frame around the heatmap."""
    for ax in cm.heatmap_axes.ravel():
        ax.grid(
            axis="both",
            which="minor",
            color="grey",
            linestyle="--",
            alpha=0.6,
            zorder=0,
        )
    if outline:
        kws = dict(outline_kws or {})
        cm.ax_heatmap.set_axis_on()
        for side in ("top", "right", "left", "bottom"):
            sp = cm.ax_heatmap.spines[side]
            sp.set_visible(True)
            sp.set_color(kws.get("color", "black"))
            sp.set_linewidth(kws.get("linewidth", 1))


# ---------------------------------------------------------------------------
# gene_dotplot2D (one gene across two grouping variables)
# ---------------------------------------------------------------------------


def _natural_sort_key(s):
    """Sort key that orders embedded numbers numerically (``E15.5 < E16 < P0 < P2 < P11``)."""
    parts = re.split(r"(\d+(?:\.\d+)?)", str(s))
    return tuple(
        (0, float(p), "") if i % 2 else (1, 0.0, p.lower()) for i, p in enumerate(parts)
    )


def _resolve_gene(var_names, gene: str) -> str:
    """Return ``gene`` as present in ``var_names``, falling back to a case-insensitive match."""
    if gene in var_names:
        return gene
    idx = pd.Index(var_names)
    hits = idx[idx.str.lower() == str(gene).lower()]
    if len(hits) == 0:
        raise KeyError(f"gene {gene!r} not found in adata.var_names.")
    print(f"gene {gene!r} not found; using case-insensitive match {hits[0]!r}.")
    return hits[0]


def gene_dotplot2D(
    adata,
    gene: str,
    row: str = "Subclass",
    col: str = "Age",
    obs=None,
    sep: str = "+",
    pseudobulk_kws: Optional[Mapping[str, Any]] = None,
    pseudobulk_adata: Optional[str] = None,
    overwrite: bool = False,
    layer: Optional[str] = None,
    min_cells: Optional[int] = 10,
    standard_scale: Optional[str] = None,
    row_order: Optional[Sequence] = None,
    col_order: Optional[Sequence] = None,
    row_cluster: bool = False,
    col_cluster: bool = False,
    cmap: str = "Reds",
    palette=None,
    legend_kws: Optional[Mapping[str, Any]] = None,
    figsize=None,
    marker: str = "o",
    plot_kws: Optional[Mapping[str, Any]] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
    outline: bool = True,
    outline_kws: Optional[Mapping[str, Any]] = None,
):
    """Dot heatmap of a single gene across two grouping variables.

    Where :func:`gene_dotplot` shows many genes across one grouping, this
    shows one gene (e.g. ``Mecp2``) across the combinations of two obs
    variables (e.g. ``Subclass`` x ``Age``): rows are the categories of
    ``row``, columns the categories of ``col``, dot colour encodes the
    expression level and dot size the fraction of expressing cells.

    Single-cell input is first aggregated with
    :func:`adataviz.tools.pseudobulk` using the composite grouping
    ``f"{row}{sep}{col}"`` (an existing ``f"{col}{sep}{row}"`` column is
    reused as well). If that column is missing from ``adata.obs`` / ``obs``
    it is built from ``row`` and ``col``. The pseudobulk AnnData is written
    to ``pseudobulk_adata`` (default ``save + ".h5ad"``); if that file
    already exists it is loaded instead and :func:`pseudobulk` is skipped
    (unless ``overwrite``), so plotting other genes is cheap.

    Parameters
    ----------
    adata : AnnData or str
        Single-cell ``.h5ad`` path (or backed AnnData), or a pseudobulk
        AnnData/path as produced by :func:`adataviz.tools.pseudobulk`
        (detected by its ``"mean"`` and ``"frac"`` layers) whose
        ``obs_names`` are ``"<row><sep><col>"`` labels.
    gene : str
        Gene to plot. Falls back to a case-insensitive match
        (``"MECP2"`` -> ``"Mecp2"``).
    row, col : str, default ``"Subclass"``, ``"Age"``
        Obs variables shown on the heatmap rows / columns. Swap them to
        transpose the plot.
    obs : DataFrame, str, or None, optional
        Alternative cell metadata (DataFrame, or TSV/CSV/``.h5ad`` path);
        None uses ``adata.obs``. Ignored for pseudobulk input.
    sep : str, default ``"+"``
        Separator joining ``row`` and ``col`` in the composite group label.
        Category names must not contain it.
    pseudobulk_kws : dict, optional
        Extra arguments for :func:`adataviz.tools.pseudobulk`. Defaults:
        ``dict(downsample=2000, use_raw=True, normalization="CPM",
        expression_cutoff=0, n_jobs=1)``.
    pseudobulk_adata : str, optional
        Path of the pseudobulk ``.h5ad`` for single-cell input. If it
        exists it is loaded and :func:`pseudobulk` is not run; otherwise
        the pseudobulk is computed and written there. Defaults to
        ``save + ".h5ad"`` (in memory only when ``save`` is also None).
        Ignored when ``adata`` is already a pseudobulk.
    overwrite : bool, default False
        Recompute the pseudobulk even if ``pseudobulk_adata`` already exists.
    layer : str, optional
        Pseudobulk matrix used for dot colour: ``"X"`` or a layer name.
        By default ``"X"`` when the pseudobulk was normalised (e.g.
        log(CPM), recorded in ``uns["Normalization"]``), else ``"mean"``.
    min_cells : int or None, default 10
        Hide groups (row/col combinations) with fewer cells than this.
    standard_scale : {None, "row", "col"}, optional
        Min-max scale the values within each row or column (e.g.
        ``"row"`` to compare ages within each subclass).
    row_order, col_order : sequence, optional
        Explicit category order (and subset) for rows / columns. Defaults
        to natural sort order (``E15.5 < E16 < P0 < P2 < P11``).
    row_cluster, col_cluster : bool, default False
        Whether to hierarchically cluster rows / columns.
    cmap : str, default ``"Reds"``
        Colormap for the dot colour.
    palette : str, dict, or None, optional
        Annotation colours: an ``.xlsx`` path with one sheet per variable
        (``Hex`` column), or a dict ``{row: {cat: color}, col: {...}}``;
        missing colours fall back to defaults.
    legend_kws : dict, optional
        Colorbar options (``cmap_legend_kws``) for ``DotClustermapPlotter``.
    figsize : tuple of float, optional
        Figure size; auto-sized from the grid when None.
    marker : str, default ``"o"``
        Dot marker.
    plot_kws : dict, optional
        Extra kwargs for :class:`PyComplexHeatmap.DotClustermapPlotter`
        (override the defaults set here).
    title : str, optional
        Figure title; defaults to the gene name.
    save : str, optional
        Figure output path; also the default pseudobulk path
        (``save + ".h5ad"``) when ``pseudobulk_adata`` is None.
    show : bool, default False
        Whether to display the figure.
    outline : bool, default True
        Draw a frame around the heatmap.
    outline_kws : dict, optional
        ``color`` / ``linewidth`` of the frame.

    Returns
    -------
    PyComplexHeatmap.DotClustermapPlotter
    """
    pch = import_pch()
    HeatmapAnnotation, anno_simple, DotClustermapPlotter = (
        pch.HeatmapAnnotation,
        pch.anno_simple,
        pch.DotClustermapPlotter,
    )
    from ..tools import load_obs, pseudobulk

    if row == col:
        raise ValueError("`row` and `col` must be different variables.")
    if standard_scale not in (None, "row", "col"):
        raise ValueError("standard_scale must be None, 'row' or 'col'.")
    plot_kws = dict(plot_kws or {})
    outline_kws = dict(outline_kws or {})

    raw = (
        anndata.read_h5ad(os.path.expanduser(adata), backed="r")
        if isinstance(adata, str)
        else adata
    )
    is_pseudobulk = {"mean", "frac"} <= set(raw.layers.keys())
    if pseudobulk_adata is not None:
        pb_path = os.path.expanduser(pseudobulk_adata)
    elif save:
        pb_path = os.path.expanduser(save) + ".h5ad"
    else:
        pb_path = None

    if is_pseudobulk:
        pb = raw
    elif pb_path is not None and os.path.exists(pb_path) and not overwrite:
        print(f"Reusing existing pseudobulk: {pb_path}")
        if isinstance(adata, str):
            raw.file.close()
        pb = anndata.read_h5ad(pb_path, backed="r")
    else:
        if isinstance(adata, str):
            adata_path = os.path.expanduser(adata)
        elif raw.isbacked:
            adata_path = str(raw.filename)
        else:
            raise ValueError(
                "Single-cell input must be an .h5ad path (or a backed AnnData): "
                "pseudobulk reads cells from disk."
            )
        cell_obs = load_obs(obs) if obs is not None else raw.obs
        if isinstance(adata, str):
            raw.file.close()
        key = next(
            (k for k in (f"{row}{sep}{col}", f"{col}{sep}{row}") if k in cell_obs.columns),
            None,
        )
        if key is None:
            missing = [c for c in (row, col) if c not in cell_obs.columns]
            if missing:
                raise KeyError(f"columns {missing} not found in obs.")
            key = f"{row}{sep}{col}"
            valid = cell_obs[row].notna() & cell_obs[col].notna()
            obs_arg = pd.DataFrame(index=cell_obs.index)
            obs_arg[key] = (
                cell_obs[row].astype(str) + sep + cell_obs[col].astype(str)
            ).where(valid)
        else:
            obs_arg = None if obs is None else cell_obs.loc[:, [key]].copy()
        kws = dict(
            downsample=2000,
            use_raw=True,
            normalization="CPM",
            expression_cutoff=0,
            n_jobs=1,
        )
        kws.update(pseudobulk_kws or {})
        print(f"Running pseudobulk with groupby={key!r} -> {pb_path or 'memory'}")
        pb = pseudobulk(adata_path, obs_path=obs_arg, groupby=key, save=pb_path, **kws)
        if pb_path is not None:
            pb = anndata.read_h5ad(pb_path, backed="r")
        if pb.obs.index.name is None:
            pb.obs.index.name = key

    gname = _resolve_gene(pb.var_names, gene)
    key_name = pb.obs.index.name or f"{row}{sep}{col}"
    sub = pb[:, [gname]]
    sub = sub.to_memory() if getattr(pb, "isbacked", False) else sub.copy()
    if getattr(pb, "isbacked", False) and pb is not adata:
        pb.file.close()

    if layer is None:
        layer = "X" if sub.uns.get("Normalization") else "mean"
    value_label = (
        str(sub.uns.get("Normalization", "Expression")) if layer == "X" else layer
    )
    if standard_scale is not None:
        value_label = f"{value_label} (scaled by {standard_scale})"
    values = sub.to_df(layer=None if layer == "X" else layer)[gname]
    frac = sub.to_df(layer="frac")[gname]

    # Split "<a><sep><b>" labels back into the two variables.
    labels = sub.obs_names.to_series()
    n_sep = labels.str.count(re.escape(sep))
    if (n_sep != 1).any():
        bad = labels[n_sep != 1].tolist()[:5]
        raise ValueError(
            f"Cannot split group labels on {sep!r} into exactly two parts, e.g. {bad}. "
            "Use a different `sep`."
        )
    key_parts = key_name.split(sep)
    if len(key_parts) == 2 and set(key_parts) == {row, col}:
        first, second = key_parts
    else:
        print(
            f"pseudobulk groupby {key_name!r} does not match row/col; "
            f"assuming labels are '<{row}>{sep}<{col}>'."
        )
        first, second = row, col
    split = labels.str.split(sep, n=1, expand=True)
    long_df = pd.DataFrame(
        {
            first: split[0].values,
            second: split[1].values,
            "Value": values.values,
            "frac": frac.values,
        },
        index=labels.index,
    )
    if min_cells and "cell_count" in sub.obs.columns:
        long_df = long_df.loc[sub.obs.loc[long_df.index, "cell_count"] >= min_cells]
    long_df = long_df.dropna(subset=["Value"])
    if long_df.empty:
        raise ValueError("No groups left to plot (check min_cells / data).")

    mat = long_df.pivot(index=row, columns=col, values="Value")
    fmat = long_df.pivot(index=row, columns=col, values="frac")
    rows = (
        [str(r) for r in row_order if str(r) in mat.index]
        if row_order is not None
        else sorted(mat.index, key=_natural_sort_key)
    )
    cols = (
        [str(c) for c in col_order if str(c) in mat.columns]
        if col_order is not None
        else sorted(mat.columns, key=_natural_sort_key)
    )
    mat, fmat = mat.loc[rows, cols], fmat.loc[rows, cols]
    if standard_scale is not None:
        ax_ = 1 if standard_scale == "row" else 0
        mn, mx = mat.min(axis=ax_), mat.max(axis=ax_)
        rng = (mx - mn).replace(0, np.nan)
        scaled = mat.sub(mn, axis=1 - ax_).div(rng, axis=1 - ax_)
        # Constant rows/columns (zero range) become 0 instead of NaN.
        mat = scaled.where(~(scaled.isna() & mat.notna()), 0.0)
    plot_data = (
        mat.rename_axis(index=row, columns=col)
        .reset_index()
        .melt(id_vars=row, var_name=col, value_name="Value")
        .merge(
            fmat.rename_axis(index=row, columns=col)
            .reset_index()
            .melt(id_vars=row, var_name=col, value_name="frac"),
            on=[row, col],
        )
        .dropna(subset=["Value"])
    )

    def _colors(var, cats) -> dict:
        pal = palette.get(var) if isinstance(palette, Mapping) else palette
        return resolve_palette(pal, cats, sheet_name=var)

    left_ha = HeatmapAnnotation(
        axis=0,
        verbose=0,
        **{
            row: anno_simple(
                pd.Series(rows, index=rows),
                colors=_colors(row, rows),
                add_text=False,
                legend=False,
                height=3,
                label=row,
            )
        },
    )
    top_ha = HeatmapAnnotation(
        axis=1,
        verbose=0,
        **{
            col: anno_simple(
                pd.Series(cols, index=cols),
                colors=_colors(col, cols),
                add_text=False,
                legend=False,
                height=3,
                label=col,
            )
        },
    )

    if legend_kws is None:
        legend_kws = dict(extendfrac=0.1, extend="both", label=value_label)
    if figsize is None:
        figsize = (max(4.0, 0.30 * len(cols) + 3.5), max(3.0, 0.24 * len(rows) + 2.0))
    if title is None:
        title = gname

    defaults = dict(
        marker=marker,
        grid=None,
        dot_legend_marker=marker,
        cmap_legend_kws=legend_kws,
        dot_legend_kws=dict(_DOT_LEGEND_KWS, title="Fraction"),
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        row_cluster_method="ward",
        row_cluster_metric="euclidean",
        col_cluster_method="ward",
        col_cluster_metric="euclidean",
        col_names_side="bottom",
        row_names_side="right",
        show_rownames=True,
        show_colnames=True,
        row_dendrogram=False,
        col_dendrogram=False,
        xticklabels_kws=dict(labelrotation=-45, labelsize=8, bottom=True),
        yticklabels_kws=dict(labelsize=8, right=True),
        legend_hpad=4,
        legend_side="right",
        legend_anchor="ax_heatmap",
        spines=False,
        # Size proportional to frac (no floor), so missing / filtered groups
        # (filled with frac=0 by PyComplexHeatmap) draw no dot.
        min_s=0,
    )
    for k, v in defaults.items():
        plot_kws.setdefault(k, v)
    if title:
        plot_kws.setdefault("xlabel", title)
        plot_kws.setdefault("xlabel_side", "top")

    fig = plt.figure(figsize=figsize)
    with quiet_pch(fig):
        cm = DotClustermapPlotter(
            data=plot_data,
            top_annotation=top_ha,
            left_annotation=left_ha,
            x_order=None if col_cluster else cols,
            y_order=None if row_cluster else rows,
            x=col,
            y=row,
            value="Value",
            c="Value",
            s="frac",
            cmap=cmap,
            verbose=0,
            **plot_kws,
        )
    _style_dot_axes(cm, outline, outline_kws)
    with quiet_pch(fig):
        save_or_show(fig, save, show=show)
    return cm


# ---------------------------------------------------------------------------
# Ridgeline (joyplot)
# ---------------------------------------------------------------------------


def _kde_curve(values, grid, bw_method=None):
    """Return a gaussian-KDE density of *values* evaluated on *grid*.

    Falls back to a normalised histogram when the values are constant (or
    too few to build a KDE), so degenerate groups still draw something.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 2 or np.allclose(v, v[0]):
        hist, edges = np.histogram(v, bins=max(10, grid.size // 8), density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        return np.interp(grid, centers, hist, left=0.0, right=0.0)
    from scipy.stats import gaussian_kde

    kde = gaussian_kde(v, bw_method=bw_method)
    return kde(grid)


def ridgeline(
    adata,
    gene: str,
    groupby: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    overlap: float = 0.7,
    bw_method: Union[None, str, float] = None,
    n_points: int = 256,
    clip: Optional[tuple] = None,
    linewidth: float = 0.8,
    linecolor: str = "white",
    fill_alpha: float = 1.0,
    label_colors: Optional[Mapping[str, str]] = None,
    figsize=None,
    ax=None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
):
    """Ridgeline plot (a.k.a. joyplot) of a gene's expression per group.

    Draws one smoothed density curve per category of ``groupby``, stacked
    vertically with a controllable amount of vertical overlap so the
    distribution of ``gene`` expression can be compared across groups at a
    glance. Densities are estimated with a gaussian KDE (falling back to a
    normalised histogram for degenerate groups) from the per-cell
    expression pulled out of ``adata`` (see ``layer`` / ``use_raw``).

    Parameters
    ----------
    adata : AnnData
        Source data. Expression for ``gene`` is read from ``adata.X``
        (or ``adata.raw.X`` / ``adata.layers[layer]``).
    gene : str
        Single gene to plot. Must exist in the resolved var space
        (``adata.raw.var_names`` when ``use_raw=True``, else
        ``adata.var_names``).
    groupby : str
        Column in ``adata.obs`` used as the categorical grouping; each
        category becomes one ridge (row).
    layer : str, optional
        Name of an ``adata.layers`` entry to read expression from. When
        given and present it takes precedence over ``adata.X``.
    use_raw : bool, default False
        If True and ``adata.raw`` exists, read expression and gene names
        from ``adata.raw``. Ignored when ``layer`` is provided.
    order : sequence, optional
        Explicit ordering of ``groupby`` categories from bottom to top.
        Categories not listed are dropped. When None the natural
        categorical order of the data is used.
    palette : dict, str, or None, optional
        Colour mapping for the ``groupby`` categories. A dict maps
        category -> colour; a string names a palette sheet resolved via
        :func:`resolve_palette`. When None, colours are resolved from
        ``adata`` (e.g. stored ``*_colors``) or a default palette.
    overlap : float, default 0.7
        Fraction of vertical overlap between adjacent ridges. ``0`` draws
        non-overlapping curves; larger values push each curve further up
        into its neighbours (the peak height equals ``1 + overlap`` rows).
    bw_method : str, float, or None, optional
        Bandwidth selector forwarded to :class:`scipy.stats.gaussian_kde`
        (e.g. ``"scott"``, ``"silverman"``, or a scalar factor). None uses
        SciPy's default.
    n_points : int, default 256
        Number of points on the shared x-grid where each density is
        evaluated.
    clip : tuple of float, optional
        ``(low, high)`` limits for the x-grid. When None the range spans
        the observed min/max of the expression (with a small margin).
    linewidth : float, default 0.8
        Width of the outline drawn on top of each filled density.
    linecolor : str, default "white"
        Colour of the density outline.
    fill_alpha : float, default 1.0
        Opacity of the filled area under each curve.
    label_colors : dict, optional
        Optional ``{category: colour}`` mapping used to colour individual
        y tick labels (e.g. to highlight one group). Categories not listed
        keep the default label colour.
    figsize : tuple of float, optional
        Figure size in inches. When None it is derived from the number of
        categories (taller for more groups).
    ax : matplotlib.axes.Axes, optional
        Existing Axes to draw into. When None a new figure and Axes are
        created.
    title : str, optional
        Axes title. When None a ``"{gene} Expression"`` title is used.
    xlabel : str, optional
        X axis label. Defaults to ``"Normalized expression"``.
    ylabel : str, optional
        Y axis label. Defaults to ``groupby``.
    save : str, optional
        Path to write the figure to (via :func:`save_or_show`).
    show : bool, default False
        Whether to display the figure (``plt.show``) via
        :func:`save_or_show`.

    Returns
    -------
    matplotlib.axes.Axes
        The Axes containing the ridgeline plot.
    """
    if not isinstance(gene, str):
        raise TypeError("ridgeline plots a single gene; pass a str.")
    long = _gene_long_df(adata, [gene], groupby, layer=layer, use_raw=use_raw)
    cats = categorical_order(long[groupby], order)
    long = long[long[groupby].isin(cats)]

    _obs, ad = resolve_adata_obs(adata)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )

    vals = long["value"].to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    if clip is None:
        lo, hi = float(np.min(vals)), float(np.max(vals))
        margin = 0.05 * (hi - lo or 1.0)
        clip = (lo - margin, hi + margin)
    grid = np.linspace(clip[0], clip[1], n_points)

    n = len(cats)
    own_fig = ax is None
    if figsize is None:
        # Reserve room for the category labels so long names don't squeeze
        # the ridges into a sliver.
        figsize = (3.2 + label_extent(cats, 8), max(3.0, 0.42 * n + 1.0))
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    step = 1.0
    amp = (1.0 + max(0.0, overlap)) * step
    densities = {}
    dmax = 0.0
    for c in cats:
        d = _kde_curve(long.loc[long[groupby] == c, "value"], grid, bw_method)
        densities[c] = d
        dmax = max(dmax, float(np.max(d)) if d.size else 0.0)
    scale = amp / dmax if dmax > 0 else 1.0

    # Draw from top to bottom so lower ridges overlap (sit in front of) the
    # ridges above them, matching the classic joyplot look.
    for i, c in enumerate(cats):
        base = i * step
        y = base + densities[c] * scale
        z = n - i  # bottom rows in front
        ax.fill_between(
            grid, base, y, color=colors[c], alpha=fill_alpha, zorder=z, lw=0
        )
        ax.plot(grid, y, color=linecolor, lw=linewidth, zorder=z + 0.5)
        ax.hlines(base, grid[0], grid[-1], color="0.3", lw=0.5, zorder=z - 0.5)

    ax.set_yticks([i * step for i in range(n)])
    ax.set_yticklabels(cats)
    if label_colors:
        for tick, c in zip(ax.get_yticklabels(), cats):
            if str(c) in {str(k) for k in label_colors}:
                tick.set_color(label_colors[c] if c in label_colors else label_colors[str(c)])
    ax.set_ylim(-0.5 * step, (n - 1) * step + amp + 0.2 * step)
    ax.set_xlim(clip[0], clip[1])
    ax.set_xlabel(xlabel if xlabel is not None else "Normalized expression")
    ax.set_ylabel(ylabel if ylabel is not None else groupby)
    ax.set_title(title if title is not None else f"{gene} Expression")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)
    apply_layout(fig, own=own_fig)
    save_or_show(fig, save, show=show)
    return ax
