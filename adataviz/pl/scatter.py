"""Embedding scatter plots: ``plot_categorical`` and ``plot_gene``.

These wrap the lower-level helpers in :mod:`adataviz.utils` (the
``categorical_scatter`` / ``continuous_scatter`` engines) so that they
can be used directly from ``adataviz.pl`` with palette resolution from
``adata.uns`` and a consistent ``legend_kws`` interface.

``plot_gene`` is now scoped to the embedding plot only; for box / strip
/ violin plots of single genes use :func:`adataviz.pl.boxplot` or
:func:`adataviz.pl.stacked_violinplot`.
"""

from __future__ import annotations

import os
from typing import Any, Mapping, Optional

import anndata
import matplotlib.pyplot as plt
import pandas as pd

from ..tools import load_adata
from ..utils import (
    categorical_scatter,
    continuous_scatter,
    normalize_mc_by_cell,
)

__all__ = ["plot_categorical", "plot_gene", "assign_default_colors"]


# ---------------------------------------------------------------------------
# Default categorical palette
# ---------------------------------------------------------------------------


_DEFAULT_20 = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
    "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
    "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
]
_DEFAULT_28 = _DEFAULT_20 + [
    "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#bd9e39",
    "#ad494a",
]


def assign_default_colors(adata, key: str) -> dict:
    """Assign default categorical colours for ``adata.obs[key]``.

    Mirrors scanpy's palette logic, writes the result to
    ``adata.uns[f"{key}_colors"]`` and returns a ``{cat: hex}`` dict.
    """
    from matplotlib.colors import to_hex
    from matplotlib import rcParams

    categories = adata.obs[key].cat.categories
    n = len(categories)
    cycle = rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(cycle) >= n:
        palette = [to_hex(c) for c in cycle[:n]]
    elif n <= 20:
        palette = _DEFAULT_20[:n]
    elif n <= 28:
        palette = _DEFAULT_28[:n]
    else:
        # Fall back to tab20 cycled
        cmap = plt.get_cmap("tab20")
        palette = [to_hex(cmap(i % cmap.N)) for i in range(n)]
    adata.uns[f"{key}_colors"] = palette
    return {cat: col for cat, col in zip(categories, palette)}


# ---------------------------------------------------------------------------
# plot_categorical
# ---------------------------------------------------------------------------


def plot_categorical(
    adata,
    ax=None,
    basis: str = "umap",
    groupby: str = "Subclass",
    coding: bool = True,
    coded_marker: bool = True,
    save: Optional[str] = None,
    palette=None,
    sheet_name: Optional[str] = None,
    show: bool = True,
    figsize=(4, 3.5),
    ncol: Optional[int] = None,
    fontsize: int = 5,
    legend_fontsize: int = 5,
    legend_kws: Optional[Mapping[str, Any]] = None,
    legend_title_fontsize: int = 5,
    marker_fontsize: int = 4,
    marker_pad: float = 0.1,
    linewidth: float = 0.5,
    axis_format: str = "tiny",
    alpha: float = 0.7,
    text_kws: Optional[Mapping[str, Any]] = None,
    **kwargs,
):
    """Categorical scatter on an embedding (e.g. UMAP) coloured by ``groupby``.

    Resolves colours from (in order) ``palette`` arg, ``adata.uns[<groupby>_colors]``,
    or a generated default palette. ``legend_kws`` is forwarded to the legend.
    """
    if basis.startswith("X_"):
        basis = basis[2:]
    if sheet_name is None:
        sheet_name = groupby
    adata = load_adata(adata)
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
        adata.obs[groupby] = adata.obs[groupby].astype("category")

    colors = None
    if palette is not None:
        if isinstance(palette, dict):
            colors = palette.copy()
        elif isinstance(palette, str):
            try:
                colors = pd.read_excel(
                    os.path.expanduser(palette), sheet_name=sheet_name, index_col=0
                ).Hex.to_dict()
                existed = adata.obs[groupby].unique().tolist()
                for k in list(colors.keys()):
                    if k not in existed:
                        del colors[k]
                for k in existed:
                    colors.setdefault(k, "gray")
            except Exception:
                colors = None
    if colors is None:
        if f"{groupby}_colors" in adata.uns:
            colors = {
                cat: col
                for cat, col in zip(
                    adata.obs[groupby].cat.categories.tolist(),
                    adata.uns[f"{groupby}_colors"],
                )
            }
        else:
            colors = assign_default_colors(adata, groupby)
    else:
        adata.uns[f"{groupby}_colors"] = [
            colors.get(k, "grey") for k in adata.obs[groupby].cat.categories.tolist()
        ]

    text_kws = {} if text_kws is None else dict(text_kws)
    text_kws.setdefault("fontsize", fontsize)
    kwargs.setdefault("hue", groupby)
    kwargs.setdefault("text_anno", groupby)
    kwargs.setdefault("text_kws", text_kws)
    kwargs.setdefault("luminance", 0.65)
    kwargs.setdefault("dodge_text", False)
    kwargs.setdefault("axis_format", axis_format)
    kwargs.setdefault("show_legend", True)
    kwargs.setdefault("marker_fontsize", marker_fontsize)
    kwargs.setdefault("marker_pad", marker_pad)
    kwargs.setdefault("linewidth", linewidth)
    kwargs.setdefault("alpha", alpha)
    kwargs["coding"] = coding
    kwargs["coded_marker"] = coded_marker

    legend_kws = {} if legend_kws is None else dict(legend_kws)
    legend_kws.setdefault("fontsize", legend_fontsize)
    legend_kws.setdefault("title", groupby)
    legend_kws.setdefault("title_fontsize", legend_title_fontsize)
    if ncol is not None:
        legend_kws.setdefault("ncol", ncol)
    kwargs.setdefault(
        "dodge_kws",
        dict(
            arrowprops=dict(
                arrowstyle="->",
                fc="grey",
                ec="none",
                connectionstyle="angle,angleA=-90,angleB=180,rad=5",
            ),
            autoalign="xy",
        ),
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=300)
    else:
        fig = ax.figure
    categorical_scatter(
        data=adata[adata.obs[groupby].notna()],
        ax=ax,
        basis=basis,
        palette=colors,
        legend_kws=legend_kws,
        **kwargs,
    )
    if save:
        fig.savefig(os.path.expanduser(save), bbox_inches="tight")
    if show:
        plt.show()
    return ax


# ---------------------------------------------------------------------------
# plot_gene (UMAP only)
# ---------------------------------------------------------------------------


def plot_gene(
    adata,
    gene: str,
    basis: str = "umap",
    obs=None,
    query_str: Optional[str] = None,
    title: Optional[str] = None,
    cmap: str = "parula",
    hue_norm=None,
    cbar_kws=None,
    axis_format: str = "tiny",
    obsm=None,
    normalize_per_cell: bool = False,
    hypo_score: bool = False,
    clip_norm_value: float = 10,
    figsize=(4, 4),
    ax=None,
    save: Optional[str] = None,
    show: bool = True,
    scatter_kws: Optional[Mapping[str, Any]] = None,
):
    """Plot a single gene's expression on an embedding.

    The previous ``plot_gene`` produced both a scatter and a per-group
    violin/box plot. The violin/box behaviour now lives in
    :func:`adataviz.pl.boxplot`; this function is purely about the
    embedding visualisation.

    Parameters
    ----------
    adata : AnnData or path
        Source data.
    gene : str
        Gene name (must exist in ``adata.var_names``).
    basis : str, default ``"umap"``
        Embedding key in ``adata.obsm`` (without ``X_`` prefix).
    obs, obsm : optional
        Override metadata or alternative embeddings (path or AnnData).
    """
    if cbar_kws is None:
        cbar_kws = dict(extendfrac=0.1)
    if scatter_kws is None:
        scatter_kws = {}

    raw = anndata.read_h5ad(os.path.expanduser(adata), backed="r") if isinstance(
        adata, str
    ) else adata
    use = raw[:, gene].to_memory() if raw.isbacked else raw[:, gene].copy()
    if isinstance(adata, str) and raw.isbacked:
        raw.file.close()

    if normalize_per_cell:
        use = normalize_mc_by_cell(
            use_adata=use,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
            hypo_score=hypo_score,
        )

    if obsm is not None:
        if isinstance(obsm, str):
            obsm = anndata.read_h5ad(os.path.expanduser(obsm), backed="r")
        keep = list(set(use.obs_names) & set(obsm.obs_names))
        use = use[keep, :]
        use.obsm = obsm[keep].obsm
        for col in obsm.obs.columns:
            if col not in use.obs.columns:
                use.obs[col] = obsm.obs.loc[use.obs_names, col].tolist()
        if hasattr(obsm, "isbacked") and obsm.isbacked:
            obsm.file.close()

    if obs is not None:
        if isinstance(obs, str):
            obs = pd.read_csv(os.path.expanduser(obs), sep="\t", index_col=0)
        else:
            obs = obs.copy()
    else:
        obs = use.obs.copy()
    if query_str is not None:
        obs = obs.query(query_str)
    overlap = list(set(use.obs_names) & set(obs.index))
    use = use[overlap, :]
    use.obs = obs.loc[use.obs_names]
    use.obs[gene] = use.to_df().loc[use.obs_names, gene].tolist()

    if title is None:
        title = query_str if query_str is not None else gene

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=300)
    else:
        fig = ax.figure
    continuous_scatter(
        data=use,
        ax=ax,
        cmap=cmap,
        hue_norm=hue_norm,
        cbar_kws=cbar_kws,
        hue=gene,
        axis_format=axis_format,
        text_anno=None,
        basis=basis if not basis.startswith("X_") else basis[2:],
        **scatter_kws,
    )
    ax.set_title(title)
    if save:
        fig.savefig(os.path.expanduser(save), bbox_inches="tight")
    if show:
        plt.show()
    return ax
