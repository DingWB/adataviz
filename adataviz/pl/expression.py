"""Per-gene expression plots that don't fit on the embedding.

- :func:`boxplot` -- multi-gene, multi-group box / violin / strip plot
  driven from a long-form pivot of ``adata.X``.
- :func:`stacked_violinplot` -- a clean seaborn-based stacked violin
  (one row per gene, one column-block per group).
- :func:`gene_dotplot` -- the rich PyComplexHeatmap dot heatmap that used
  to live as ``plot_genes`` (supports ``parent_col``, ``query_str``, etc.).
- :func:`get_genes_mean_frac` -- helper used by both ``gene_dotplot`` and
  the interactive variant.
"""

from __future__ import annotations

import os
from typing import Any, Mapping, Optional, Sequence, Union

import anndata
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
from ..utils import normalize_mc_by_cell

__all__ = [
    "boxplot",
    "stacked_violinplot",
    "gene_dotplot",
    "get_genes_mean_frac",
]


# ---------------------------------------------------------------------------
# Multi-gene / multi-group boxplot (or violin / strip)
# ---------------------------------------------------------------------------


def _gene_long_df(adata, genes, groupby, layer=None, use_raw=False):
    """Return a tidy ``[groupby, gene, value]`` DataFrame."""
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("boxplot/stacked_violinplot require an AnnData.")
    src = adata.raw if use_raw and adata.raw is not None else adata
    var = list(src.var_names)
    missing = [g for g in genes if g not in var]
    if missing:
        raise KeyError(f"genes not in adata.var_names: {missing[:5]}")
    idx = [var.index(g) for g in genes]
    if layer is not None and layer in adata.layers:
        X = adata.layers[layer][:, idx]
    else:
        X = src.X[:, idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    df = pd.DataFrame(np.asarray(X), index=adata.obs_names, columns=list(genes))
    df[groupby] = adata.obs[groupby].astype(str).values
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

    Parameters
    ----------
    adata : AnnData
        Source data.
    genes : str or list of str
        Gene(s) to plot. A single gene draws one panel; multiple genes
        draw one panel per gene (sharing the X axis).
    groupby : str
        Column in ``adata.obs`` used as the X axis.
    kind : {"violin", "box", "strip"}, default ``"violin"``
        Plot style. ``"violin"`` uses ``seaborn.violinplot`` with
        ``cut=0`` (so distributions don't bleed past the data range).
    show_strip : bool, default False
        Overlay a low-alpha strip plot on top of the violin/box.
    legend_kws : dict, optional
        Forwarded to :meth:`matplotlib.axes.Axes.legend`.
    """
    import seaborn as sns

    if isinstance(genes, str):
        genes = [genes]
    long = _gene_long_df(adata, genes, groupby, layer=layer, use_raw=use_raw)
    cats = categorical_order(long[groupby], order)
    long = long[long[groupby].isin(cats)]
    long[groupby] = pd.Categorical(long[groupby], categories=cats, ordered=True)

    _obs, ad = resolve_adata_obs(adata)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
    pal = [colors[c] for c in cats]

    if gene_order is None:
        gene_order = list(genes)

    n_genes = len(gene_order)
    if figsize is None:
        figsize = (max(5, len(cats) * 0.5), 2.5 * n_genes if n_genes > 1 else 3.5)
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

    for i, g in enumerate(gene_order):
        sub = long[long["gene"] == g]
        a = axes[i]
        if kind == "violin":
            sns.violinplot(
                data=sub,
                x=groupby,
                y="value",
                order=cats,
                palette=pal,
                ax=a,
                cut=0,
                density_norm="width",
                linewidth=0.6,
                inner="quart",
            )
        elif kind == "box":
            sns.boxplot(
                data=sub,
                x=groupby,
                y="value",
                order=cats,
                palette=pal,
                ax=a,
                showfliers=False,
                linewidth=0.6,
                saturation=0.85,
            )
        elif kind == "strip":
            sns.stripplot(
                data=sub,
                x=groupby,
                y="value",
                order=cats,
                palette=pal,
                ax=a,
                size=strip_size,
                edgecolor="white",
                linewidth=0,
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
        if i < n_genes - 1:
            a.set_xlabel("")
            a.tick_params(labelbottom=False)
        else:
            a.set_xlabel(groupby)
            plt.setp(
                a.get_xticklabels(),
                rotation=rotation,
                ha="left" if rotation < 0 else "right",
            )
        if show_legend and i == 0 and legend_kws is not None:
            from matplotlib.patches import Patch

            handles = [Patch(color=colors[c], label=c) for c in cats]
            a.legend(handles=handles, **merge_legend_kws(legend_kws, title=groupby))
    if title:
        fig.suptitle(title)
    fig.tight_layout()
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
    rotation: float = -45,
    inner: Optional[str] = "quart",
):
    """Compact stacked violin plot.

    Each row (or column with ``swap_axes=True``) is one gene; each x-tick
    is a group. ``standardize='gene'`` rescales each gene to ``[0, 1]`` so
    relative differences across groups are comparable across genes.

    No scanpy dependency.
    """
    import seaborn as sns

    long = _gene_long_df(adata, genes, groupby, layer=layer, use_raw=use_raw)
    cats = categorical_order(long[groupby], order)
    long = long[long[groupby].isin(cats)]
    long[groupby] = pd.Categorical(long[groupby], categories=cats, ordered=True)

    if standardize == "gene":

        def _norm(v):
            mn, mx = v.min(), v.max()
            return (v - mn) / (mx - mn) if mx > mn else v * 0

        long["value"] = long.groupby("gene", observed=True)["value"].transform(_norm)
    elif standardize == "group":

        def _norm(v):
            mn, mx = v.min(), v.max()
            return (v - mn) / (mx - mn) if mx > mn else v * 0

        long["value"] = long.groupby([groupby, "gene"], observed=True)[
            "value"
        ].transform(_norm)

    _obs, ad = resolve_adata_obs(adata)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
    pal = [colors[c] for c in cats]

    n_genes = len(genes)
    if swap_axes:
        if figsize is None:
            figsize = (max(4, len(cats) * 0.4), max(3, n_genes * 0.4))
        fig, axes = plt.subplots(1, n_genes, figsize=figsize, sharey=True)
        axes = list(np.atleast_1d(axes).ravel())
        for i, g in enumerate(genes):
            sub = long[long["gene"] == g]
            sns.violinplot(
                data=sub,
                y=groupby,
                x="value",
                order=cats,
                palette=pal,
                ax=axes[i],
                cut=0,
                density_norm="width",
                linewidth=0.5,
                inner=inner,
                orient="h",
            )
            axes[i].set_xlabel(g, rotation=rotation, ha="right")
            axes[i].set_ylabel("")
            axes[i].set_xticks([])
            axes[i].spines["bottom"].set_visible(False)
    else:
        if figsize is None:
            figsize = (max(4, len(cats) * 0.4), max(2.5, n_genes * 0.6))
        fig, axes = plt.subplots(n_genes, 1, figsize=figsize, sharex=True)
        axes = list(np.atleast_1d(axes).ravel())
        for i, g in enumerate(genes):
            sub = long[long["gene"] == g]
            sns.violinplot(
                data=sub,
                x=groupby,
                y="value",
                order=cats,
                palette=pal,
                ax=axes[i],
                cut=0,
                density_norm="width",
                linewidth=0.5,
                inner=inner,
            )
            axes[i].set_ylabel(g, rotation=0, ha="right", va="center", labelpad=15)
            axes[i].set_yticks([])
            axes[i].spines["left"].set_visible(False)
            if i < n_genes - 1:
                axes[i].set_xlabel("")
                axes[i].tick_params(labelbottom=False)
            else:
                axes[i].set_xlabel(groupby)
                plt.setp(
                    axes[i].get_xticklabels(),
                    rotation=rotation,
                    ha="left" if rotation < 0 else "right",
                )
    if title:
        fig.suptitle(title)
    fig.tight_layout()
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

    Returns a tidy DataFrame ``[groupby, "Gene", "Mean", "frac"]``.

    Parameters
    ----------
    query : str, optional
        Pandas-style query string applied to ``adata.obs`` to subset cells
        before loading expression data into memory. Useful with backed
        AnnData to avoid materialising the full matrix.
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
        # doesn't carry one). Use ``stack(dropna=False)`` so genes whose
        # mean is NaN for a queried single group don't silently drop the
        # entire row.
        mean_df = use.to_df(layer=layer)
        frac_df = use.to_df(layer="frac")

        plot_data = mean_df.stack(dropna=False).reset_index()
        plot_data.columns = [groupby, "Gene", "Mean"]

        frac_long = frac_df.stack(dropna=False).reset_index()
        frac_long.columns = [groupby, "Gene", "frac"]

        plot_data = plot_data.merge(frac_long, on=[groupby, "Gene"], how="left")
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
    overlap = list(set(use.obs_names) & set(obs.index))
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
    if modality in ("RNA", "ATAC") and isinstance(expression_cutoff, str):
        if expression_cutoff == "median":
            cutoff = df.stack().median()
        elif expression_cutoff == "mean":
            cutoff = df.stack().mean()
        else:
            cutoff = df.stack().quantile(
                float(expression_cutoff.replace("p", "")) / 100
            )
        expression_cutoff = cutoff
    df[groupby] = obs.loc[df.index, groupby].tolist()
    plot_data = df.groupby(groupby, observed=True).mean().stack().reset_index()
    plot_data.columns = [groupby, "Gene", "Mean"]
    if "frac" in use.layers:
        D = use.to_df(layer="frac").stack().to_dict()
    elif modality in ("RNA", "ATAC"):
        frac = df.groupby(groupby, observed=True).agg(
            lambda x: (x > expression_cutoff).sum() / len(x)
        )
        D = frac.stack().to_dict()
    else:
        hypo = df.groupby(groupby, observed=True).agg(lambda x: (x < 1).sum() / len(x))
        D = hypo.stack().to_dict()
    plot_data["frac"] = (
        plot_data.loc[:, [groupby, "Gene"]]
        .apply(lambda x: tuple(x.tolist()), axis=1)
        .map(D)
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
    figsize=(10, 4),
    marker: str = "o",
    plot_kws: Optional[Mapping[str, Any]] = None,
    transpose: bool = False,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
):
    """Rich PyComplexHeatmap dot heatmap with hierarchical annotations.

    Compared to :func:`adataviz.pl.complex_dotplot`, this function adds
    optional ``parent_col`` (for two-level annotation), ``query_str``
    pre-filtering, and full ``plot_kws`` passthrough to
    ``DotClustermapPlotter``. Useful for figure-quality marker panels.

    All extra ``plot_kws`` are forwarded directly to
    :class:`PyComplexHeatmap.DotClustermapPlotter`.
    """
    from PyComplexHeatmap import (
        HeatmapAnnotation,
        anno_simple,
        DotClustermapPlotter,
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
    use = raw[:, keep].to_memory()
    if use_raw and use.raw is not None:
        rraw = use.raw[:, use.var_names].to_adata()
        use.X = rraw[use.obs_names, use.var_names].X.copy()
    if isinstance(adata, str) and raw.isbacked:
        raw.file.close()

    obs = _resolve_obs_arg(use, obs)
    if query_str is not None:
        obs = obs.query(query_str)
    overlap = list(set(use.obs_names) & set(obs.index))
    obs = obs.loc[overlap]
    use = use[overlap, :]

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

    df = use.to_df()
    if modality in ("RNA", "ATAC") and isinstance(expression_cutoff, str):
        if expression_cutoff == "median":
            cutoff = df.stack().median()
        elif expression_cutoff == "mean":
            cutoff = df.stack().mean()
        else:
            cutoff = df.stack().quantile(
                float(expression_cutoff.replace("p", "")) / 100
            )
        expression_cutoff = cutoff
    df[groupby] = use.obs.loc[df.index, groupby].tolist()
    if parent_col is not None:
        group2parent = (
            use.obs.loc[:, [groupby, parent_col]]
            .drop_duplicates()
            .set_index(groupby)[parent_col]
            .to_dict()
        )
    plot_data = df.groupby(groupby, observed=True).mean().stack().reset_index()
    plot_data.columns = [groupby, "Gene", "Mean"]
    if "frac" in use.layers:
        D = use.to_df(layer="frac").stack().to_dict()
    elif modality in ("RNA", "ATAC"):
        frac = df.groupby(groupby, observed=True).agg(
            lambda x: (x > expression_cutoff).sum() / len(x)
        )
        D = frac.stack().to_dict()
    else:
        hypo = df.groupby(groupby, observed=True).agg(lambda x: (x < 1).sum() / len(x))
        D = hypo.stack().to_dict()
    plot_data["frac"] = (
        plot_data.loc[:, [groupby, "Gene"]]
        .apply(lambda x: tuple(x.tolist()), axis=1)
        .map(D)
    )

    df_cols = pd.DataFrame(sorted(use.obs[groupby].unique()), columns=[groupby])
    if parent_col is not None:
        df_cols[parent_col] = df_cols[groupby].map(group2parent)
        df_cols.sort_values([parent_col, groupby], inplace=True)
    df_cols.index = df_cols[groupby].tolist()
    if cell_type_order is not None:
        df_cols = df_cols.loc[[c for c in cell_type_order if c in df_cols.index]]
    df_cols.dropna(inplace=True)

    axis = 1 if not transpose else 0

    if "+" in groupby:
        anno_dict = {}
        for sub in groupby.split("+"):
            df_cols[sub] = df_cols[groupby].apply(
                lambda x, _s=sub: x.split("+")[groupby.split("+").index(_s)]
            )
            sub_colors = {k: color_palette[sub][k] for k in df_cols[sub].unique()}
            anno_dict[sub] = anno_simple(
                df_cols[sub],
                colors=sub_colors,
                add_text=False,
                legend=False,
                height=3,
                label=sub,
            )
        if parent_col is not None:
            pcolors = {
                k: color_palette[parent_col][k] for k in df_cols[parent_col].unique()
            }
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
        gcolors = {k: color_palette[groupby][k] for k in df_cols[groupby].unique()}
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
            pcolors = {
                k: color_palette[parent_col][k] for k in df_cols[parent_col].unique()
            }
            kw["parent"] = anno_simple(
                df_cols[parent_col],
                colors=pcolors,
                add_text=False,
                legend=True,
                height=3,
                label=parent_col,
            )
        col_ha = HeatmapAnnotation(**kw)

    # Always show tick labels; pick a side that does not collide with col_ha.
    # When col_ha is on top (axis=1, default) we keep group names at the
    # bottom of the heatmap and gene names on the right - both sides of the
    # heatmap stay clean and labels are clearly readable.
    if not transpose:
        top_anno, left_anno = col_ha, None
        x, y = groupby, "Gene"
        x_order = df_cols.index.tolist()
        y_order = gene_order
        col_side, row_side = "bottom", "left"
    else:
        top_anno, left_anno = None, col_ha
        x, y = "Gene", groupby
        x_order = gene_order
        y_order = df_cols.index.tolist()
        col_side, row_side = "bottom", "right"

    defaults = dict(
        marker=marker,
        grid=None,
        dot_legend_marker=marker,
        cmap_legend_kws=legend_kws,
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

    fig = plt.figure(figsize=figsize)
    cm = DotClustermapPlotter(
        data=plot_data,
        top_annotation=top_anno,
        left_annotation=left_anno,
        x_order=x_order,
        y_order=y_order,
        x=x,
        y=y,
        value="Mean",
        c="Mean",
        s="frac",
        cmap=cmap,
        verbose=0,
        **plot_kws,
    )
    for ax in cm.heatmap_axes.ravel():
        ax.grid(
            axis="both",
            which="minor",
            color="grey",
            linestyle="--",
            alpha=0.6,
            zorder=0,
        )
    _strip_cbar_white_lines(cm)
    if title:
        fig.suptitle(title)
    save_or_show(fig, save, show=show)
    return cm
