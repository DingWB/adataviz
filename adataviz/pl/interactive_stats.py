"""Interactive (Plotly) versions of the static cell-statistics plots.

These functions return :class:`plotly.graph_objects.Figure` instances
that can be displayed inline in Jupyter or saved to HTML. They mirror
the static API in :mod:`adataviz.pl.composition` and
:mod:`adataviz.pl.flow`.

All functions lazy-import :mod:`plotly` so that ``import adataviz`` does
not require it.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd

from ._utils import (
    categorical_order,
    resolve_adata_obs,
    resolve_palette,
)

__all__ = [
    "interactive_rose",
    "interactive_ring",
    "interactive_pie",
    "interactive_area",
    "interactive_trend",
    "interactive_sankey",
]


def _maybe_save(fig, save: Optional[str]):
    """Save a Plotly figure as HTML or PNG depending on extension."""
    if not save:
        return
    import os

    out = os.path.abspath(os.path.expanduser(save))
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    if out.endswith(".html"):
        fig.write_html(out, include_plotlyjs="cdn")
    else:
        fig.write_image(out)


# ---------------------------------------------------------------------------
# Rose
# ---------------------------------------------------------------------------


def interactive_rose(
    adata,
    groupby: str,
    order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
):
    """Plotly polar bar chart (Nightingale rose) of category counts."""
    import plotly.graph_objects as go

    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
    fig = go.Figure(
        go.Barpolar(
            r=counts.values,
            theta=cats,
            marker_color=[colors[c] for c in cats],
            marker_line_color="white",
            marker_line_width=1,
            opacity=0.85,
            hovertemplate="%{theta}: %{r}<extra></extra>",
        )
    )
    fig.update_layout(
        title=title,
        polar=dict(
            radialaxis=dict(showticklabels=False, ticks=""),
            angularaxis=dict(direction="clockwise"),
        ),
        height=height,
        margin=dict(l=20, r=20, t=40, b=20),
    )
    _maybe_save(fig, save)
    return fig


# ---------------------------------------------------------------------------
# Ring (donut)
# ---------------------------------------------------------------------------


def interactive_ring(
    adata,
    groupby: str,
    order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    hole: float = 0.45,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
):
    """Plotly donut chart of ``groupby`` proportions."""
    import plotly.graph_objects as go

    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
    fig = go.Figure(
        go.Pie(
            labels=cats,
            values=counts.values,
            hole=hole,
            marker=dict(colors=[colors[c] for c in cats]),
            textinfo="percent",
            hovertemplate="%{label}: %{value} (%{percent})<extra></extra>",
        )
    )
    fig.update_layout(title=title, height=height, margin=dict(l=20, r=20, t=40, b=20))
    _maybe_save(fig, save)
    return fig


# ---------------------------------------------------------------------------
# Pie
# ---------------------------------------------------------------------------


def interactive_pie(
    adata,
    groupby: str,
    order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
):
    """Plotly pie chart of ``groupby`` proportions (no hole)."""
    return interactive_ring(
        adata,
        groupby=groupby,
        order=order,
        palette=palette,
        hole=0.0,
        title=title,
        save=save,
        height=height,
    )


# ---------------------------------------------------------------------------
# Stacked area
# ---------------------------------------------------------------------------


def interactive_area(
    adata,
    groupby: str,
    split_by: str,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    normalize: bool = True,
    palette: Union[None, Mapping[str, str], str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
):
    """Plotly stacked area chart, mirroring :func:`area_plot`."""
    import plotly.graph_objects as go

    obs, ad = resolve_adata_obs(adata)
    splits = categorical_order(obs[split_by], split_order)
    cats = categorical_order(obs[groupby], order)
    ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
    ct = ct.reindex(index=cats, columns=splits, fill_value=0).astype(float)
    if normalize:
        col_sums = ct.sum(axis=0).replace(0, np.nan)
        ct = ct.div(col_sums, axis=1).fillna(0)

    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)
    fig = go.Figure()
    for c in cats:
        fig.add_trace(
            go.Scatter(
                x=splits,
                y=ct.loc[c].values,
                mode="lines",
                stackgroup="one",
                name=c,
                line=dict(width=0.5, color=colors[c]),
                fillcolor=colors[c],
                hovertemplate=f"{c}: %{{y:.3f}}<extra></extra>",
            )
        )
    fig.update_layout(
        title=title,
        height=height,
        xaxis_title=split_by,
        yaxis_title="Fraction" if normalize else "Count",
        margin=dict(l=40, r=20, t=40, b=40),
    )
    _maybe_save(fig, save)
    return fig


# ---------------------------------------------------------------------------
# Trend (line plot)
# ---------------------------------------------------------------------------


def interactive_trend(
    adata,
    groupby: str,
    split_by: str,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    normalize: bool = True,
    palette: Union[None, Mapping[str, str], str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
):
    """Plotly line trend mirror of :func:`trend_plot`."""
    import plotly.graph_objects as go

    obs, ad = resolve_adata_obs(adata)
    splits = categorical_order(obs[split_by], split_order)
    cats = categorical_order(obs[groupby], order)
    ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
    ct = ct.reindex(index=cats, columns=splits, fill_value=0).astype(float)
    if normalize:
        col_sums = ct.sum(axis=0).replace(0, np.nan)
        ct = ct.div(col_sums, axis=1).fillna(0)
    colors = resolve_palette(palette, cats, sheet_name=groupby, adata=ad, groupby=groupby)

    fig = go.Figure()
    for c in cats:
        fig.add_trace(
            go.Scatter(
                x=splits,
                y=ct.loc[c].values,
                mode="lines+markers",
                name=c,
                line=dict(color=colors[c]),
                marker=dict(size=6),
            )
        )
    fig.update_layout(
        title=title,
        height=height,
        xaxis_title=split_by,
        yaxis_title="Fraction" if normalize else "Count",
        margin=dict(l=40, r=20, t=40, b=40),
    )
    _maybe_save(fig, save)
    return fig


# ---------------------------------------------------------------------------
# Sankey
# ---------------------------------------------------------------------------


def interactive_sankey(
    adata,
    left: str,
    right: str,
    left_order: Optional[Sequence] = None,
    right_order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
):
    """Plotly Sankey diagram mirroring :func:`sankey_plot`.

    Each link can be hovered for the underlying cell count.
    """
    import plotly.graph_objects as go

    obs, ad = resolve_adata_obs(adata)
    left_cats = categorical_order(obs[left], left_order)
    right_cats = categorical_order(obs[right], right_order)
    ct = pd.crosstab(obs[left].astype(str), obs[right].astype(str))
    ct = ct.reindex(index=left_cats, columns=right_cats, fill_value=0)

    # Drop empty nodes (mirrors the static fix).
    left_cats = [c for c in left_cats if ct.loc[c].sum() > 0]
    right_cats = [c for c in right_cats if ct[c].sum() > 0]
    ct = ct.loc[left_cats, right_cats]

    nodes = list(left_cats) + list(right_cats)
    node_index = {n: i for i, n in enumerate(nodes[: len(left_cats)])}
    for i, n in enumerate(right_cats):
        node_index[(right, n)] = len(left_cats) + i

    sources, targets, values = [], [], []
    for i, lc in enumerate(left_cats):
        for j, rc in enumerate(right_cats):
            v = ct.iat[i, j]
            if v <= 0:
                continue
            sources.append(node_index[lc])
            targets.append(node_index[(right, rc)])
            values.append(int(v))

    left_colors = resolve_palette(palette, left_cats, sheet_name=left, adata=ad, groupby=left)
    right_colors = resolve_palette(None, right_cats, sheet_name=right, adata=ad, groupby=right)
    node_colors = [left_colors[c] for c in left_cats] + [
        right_colors[c] for c in right_cats
    ]

    def _hex_to_rgba(hex_color: str, alpha: float = 0.4) -> str:
        h = hex_color.lstrip("#")
        r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
        return f"rgba({r},{g},{b},{alpha})"

    link_colors = [_hex_to_rgba(left_colors[left_cats[s]]) for s in sources]

    fig = go.Figure(
        go.Sankey(
            arrangement="snap",
            node=dict(
                label=nodes,
                pad=12,
                thickness=14,
                color=node_colors,
                line=dict(color="white", width=0.5),
            ),
            link=dict(
                source=sources,
                target=targets,
                value=values,
                color=link_colors,
            ),
        )
    )
    fig.update_layout(
        title=title, height=height, margin=dict(l=20, r=20, t=40, b=20)
    )
    _maybe_save(fig, save)
    return fig


# ---------------------------------------------------------------------------
# Embedding (categorical or continuous gene) -- migrated from legacy
# ---------------------------------------------------------------------------


def interactive_embedding(
    adata=None,
    obs=None,
    variable=None,
    gene=None,
    annotation=None,
    basis="umap",
    vmin="p1",
    vmax="p99",
    cmap="jet",
    title=None,
    width=900,
    height=750,
    palette=None,
    size=None,
    show=True,
    downsample=None,
    target_fill=0.05,
    normalize_per_cell=False,
    clip_norm_value=10,
    renderer="notebook",
):
    """Interactive Plotly scatter on an embedding (UMAP/etc.).

    Pass ``variable`` for a categorical colouring or ``gene`` for a
    continuous one. ``annotation`` adds an extra hover field.
    """
    import plotly.express as px
    import plotly.io as pio

    from ._utils import show_fig
    from ..tools import load_adata, load_obs, load_color_palette
    from ..utils import normalize_mc_by_cell

    if renderer is not None:
        pio.renderers.default = renderer
    use_col = gene if gene is not None else variable
    if gene is not None and adata is None:
        raise ValueError("`gene` provided, `adata` must be provided too.")
    if adata is not None:
        adata = load_adata(adata)
    use_adata = None
    if gene is not None:
        if adata.isbacked:
            use_adata = adata[:, gene].to_memory()
        else:
            use_adata = adata[:, gene].copy()
        if normalize_per_cell:
            use_adata = normalize_mc_by_cell(
                use_adata=use_adata, normalize_per_cell=normalize_per_cell,
                clip_norm_value=clip_norm_value, hypo_score=False,
            )
    elif adata is not None:
        use_adata = adata
    if obs is None and use_adata is None:
        raise ValueError("Either `adata` or `obs` must be provided.")
    if obs is None:
        obs = use_adata.obs.copy()
    else:
        obs = load_obs(obs)
        if use_adata is not None:
            overlap = obs.index.intersection(use_adata.obs_names)
            obs = obs.loc[overlap]
            use_adata = use_adata[overlap, :]
    if gene is not None:
        obs[gene] = use_adata.to_df()[gene].tolist()
    if f"X_{basis}" in use_adata.obsm:
        obs[f"{basis}_0"] = np.asarray(use_adata.obsm[f"X_{basis}"][:, 0])
        obs[f"{basis}_1"] = np.asarray(use_adata.obsm[f"X_{basis}"][:, 1])
    if adata is not None and adata.isbacked:
        adata.file.close()

    n_points = obs.shape[0]
    if downsample is not None and n_points > downsample:
        sample_idx = np.random.choice(n_points, size=downsample, replace=False)
        obs = obs.iloc[sample_idx]

    if obs.dtypes[use_col] not in ["object", "category"]:
        vmin_q = float(int(vmin.replace("p", "")) / 100)
        vmax_q = float(int(vmax.replace("p", "")) / 100)
        range_color = [obs[use_col].quantile(vmin_q), obs[use_col].quantile(vmax_q)]
        color_discrete_map = None
        category_orders = None
    else:
        if isinstance(palette, str) and palette.endswith(".xlsx"):
            color_discrete_map = load_color_palette(
                palette=palette, adata=use_adata, groups=use_col
            )
        else:
            color_discrete_map = palette
        if color_discrete_map is not None:
            keys = list(color_discrete_map.keys())
            for k in keys:
                if k not in obs[use_col].unique().tolist():
                    del color_discrete_map[k]
        range_color = None
        order = sorted(obs[use_col].dropna().unique().tolist())
        obs[use_col] = pd.Categorical(obs[use_col], categories=order, ordered=True)
        category_orders = {use_col: order}

    keep_cols = ["cell", f"{basis}_0", f"{basis}_1"]
    if variable is not None:
        keep_cols.append(variable)
    if gene is not None:
        keep_cols.append(gene)
    if annotation is not None:
        keep_cols.append(annotation)
    obs = obs.reset_index(names="cell").loc[:, keep_cols]

    hover_data = {"cell": True, f"{basis}_0": ":0.3f", f"{basis}_1": ":0.3f"}
    if variable is not None:
        hover_data[variable] = True
    if annotation is not None:
        hover_data[annotation] = True
    if gene is not None:
        hover_data[gene] = ":.3f"
    fig = px.scatter(
        obs, x=f"{basis}_0", y=f"{basis}_1", color=use_col,
        category_orders=category_orders, hover_data=hover_data,
        range_color=range_color,
        color_discrete_sequence=px.colors.qualitative.D3,
        color_discrete_map=color_discrete_map,
        color_continuous_scale=cmap,
        template="plotly_white", render_mode="webgl",
    )
    fig.update_xaxes(
        range=[obs[f"{basis}_0"].min() - 0.5, obs[f"{basis}_0"].max() + 0.5],
        tickfont_size=12,
    )
    fig.update_yaxes(
        range=[obs[f"{basis}_1"].min() - 0.5, obs[f"{basis}_1"].max() + 0.5],
        tickfont_size=12,
    )
    if size is None:
        marker_diam_area = 2 * np.sqrt((width * height * target_fill) / (np.pi * n_points))
        marker_diam_log = 16 - 2 * np.log10(n_points)
        marker_diam = 0.7 * marker_diam_area + 0.5 * marker_diam_log
        size = int(np.clip(marker_diam, 4, 20))
    opacity = 0.8 if n_points < 500_000 else 0.6
    fig.update_traces(
        marker=dict(size=size, opacity=opacity, line=dict(width=0.12, color="black")),
        selector=dict(mode="markers"),
    )
    if title is None:
        title = f"{basis.upper()} (Colored by {use_col})"
    fig.update_layout(
        title=dict(text=title, font_size=16, x=0.5, pad=dict(t=10)),
        xaxis_title=f"{basis}_0".upper(),
        yaxis_title=f"{basis}_1".upper(),
        autosize=True, width=width, height=height, legend_title=use_col,
        legend=dict(font_size=12, itemsizing="constant", itemwidth=30, borderwidth=0.1),
    )
    if show:
        show_fig(fig, filename=f"{basis}.{use_col}")
        return None
    return fig


def interactive_categorical(adata, groupby, basis="umap", **kwargs):
    """Convenience wrapper: interactive scatter coloured by a categorical column."""
    return interactive_embedding(adata=adata, variable=groupby, basis=basis, **kwargs)


def interactive_gene(adata, gene, basis="umap", **kwargs):
    """Convenience wrapper: interactive scatter coloured by a gene's expression."""
    return interactive_embedding(adata=adata, gene=gene, basis=basis, **kwargs)


# ---------------------------------------------------------------------------
# Dot heatmap
# ---------------------------------------------------------------------------


def interactive_dotHeatmap(
    adata=None,
    obs=None,
    genes=None,
    groupby="Subclass",
    modality="RNA",
    title=None,
    use_raw=False,
    expression_cutoff="p5",
    normalize_per_cell=True,
    clip_norm_value=10,
    width=900,
    height=700,
    gene_order=None,
    colorscale="greens",
    vmin="p1",
    vmax="p99",
    show=True,
    reversescale=False,
    size_min=5,
    size_max=30,
    renderer="notebook",
):
    """Interactive Plotly dot heatmap: dot size = frac, colour = mean."""
    import plotly.io as pio
    import plotly.graph_objects as go

    from ._utils import show_fig
    from .expression import get_genes_mean_frac

    if renderer is not None:
        pio.renderers.default = renderer
    plot_data = get_genes_mean_frac(
        adata, obs=obs, groupby=groupby, modality=modality, use_raw=use_raw,
        expression_cutoff=expression_cutoff, genes=genes,
        normalize_per_cell=normalize_per_cell, clip_norm_value=clip_norm_value,
        hypo_score=False,
    )
    x_labels = plot_data[groupby].unique().tolist()
    if gene_order is None:
        y_labels = plot_data["Gene"].unique().tolist()
    else:
        y_labels = [g for g in gene_order if g in plot_data["Gene"].unique()]
    plot_data["x_cat"] = pd.Categorical(plot_data[groupby], categories=x_labels)
    plot_data["y_cat"] = pd.Categorical(plot_data["Gene"], categories=y_labels)
    frac_vals = plot_data["frac"].fillna(0).astype(float)
    sizes = (frac_vals * (size_max - size_min) + size_min).tolist()
    mean_vals = plot_data["Mean"].astype(float).tolist()
    hover_text = [
        f"Group: {g}<br>Gene: {ge}<br>Mean: {m:.4g}<br>Frac: {f:.3g}"
        for g, ge, m, f in zip(
            plot_data[groupby].tolist(), plot_data["Gene"].tolist(),
            mean_vals, frac_vals,
        )
    ]
    vmin_q = float(int(vmin.replace("p", "")) / 100)
    vmax_q = float(int(vmax.replace("p", "")) / 100)
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=plot_data[groupby].tolist(), y=plot_data["Gene"].tolist(),
        mode="markers",
        marker=dict(
            size=sizes, color=mean_vals, colorscale=colorscale,
            showscale=True, colorbar=dict(title="Mean"),
            reversescale=reversescale, sizemode="area", opacity=0.9,
            cmin=plot_data["Mean"].quantile(vmin_q),
            cmax=plot_data["Mean"].quantile(vmax_q),
        ),
        text=hover_text, hoverinfo="text",
    ))
    fig.update_xaxes(type="category", categoryorder="array",
                     categoryarray=x_labels, tickangle=-45)
    fig.update_yaxes(type="category", categoryorder="array",
                     categoryarray=list(reversed(y_labels)))
    if title is None:
        title = groupby
    fig.update_layout(
        title=title or "", xaxis_title=groupby, yaxis_title="Gene",
        width=width, height=height, plot_bgcolor="white",
    )
    if show:
        show_fig(fig, filename=f"dotHeatmap.{groupby}")
        return None
    return fig


# ---------------------------------------------------------------------------
# Boxplot
# ---------------------------------------------------------------------------


def _has_stats(adata):
    import anndata
    if isinstance(adata, str):
        adata = anndata.read_h5ad(adata, backed="r")
    return all(k in adata.layers for k in ["min", "q25", "q50", "q75", "max", "mean", "std"])


def _get_boxplot_data(adata, variable, gene, obs=None):
    import os
    import anndata

    assert isinstance(adata, anndata.AnnData)
    if adata.isbacked:
        use_adata = adata[:, gene].to_memory()
    else:
        use_adata = adata[:, gene].copy()
    if isinstance(obs, str):
        obs_path = os.path.abspath(os.path.expanduser(obs))
        sep = "\t" if obs_path.endswith((".tsv", ".txt")) else ","
        data = pd.read_csv(obs_path, index_col=0, sep=sep)
    else:
        assert isinstance(obs, pd.DataFrame)
        data = obs.copy()
    overlap = data.index.intersection(use_adata.obs_names)
    data = data.loc[overlap]
    use_adata = use_adata[overlap, :]
    data[gene] = use_adata.to_df()[gene].tolist()
    return data.loc[:, [variable, gene]]


def interactive_boxplot(
    adata,
    variable,
    gene,
    obs=None,
    palette=None,
    title=None,
    width=1100,
    height=700,
    show=True,
    renderer="notebook",
):
    """Interactive Plotly boxplot of a gene by a categorical variable.

    Auto-detects whether ``adata`` carries pre-computed stat layers
    (pseudobulk). For multi-gene support use :func:`interactive_violin`.
    """
    import anndata
    import plotly.express as px
    import plotly.io as pio
    import plotly.graph_objects as go

    from ._utils import show_fig
    from ..tools import load_color_palette

    if renderer is not None:
        pio.renderers.default = renderer
    if isinstance(adata, str):
        adata = anndata.read_h5ad(adata, backed="r")
    if obs is None:
        obs = adata.obs.copy()
    if not _has_stats(adata):
        plot_df = _get_boxplot_data(adata, variable, gene, obs=obs)
        range_y = [plot_df[gene].quantile(0.01), plot_df[gene].quantile(0.99)]
        cmap = load_color_palette(palette=palette, adata=adata, groups=variable)
        if cmap is not None:
            for k in list(cmap.keys()):
                if k not in plot_df[variable].unique().tolist():
                    del cmap[k]
        if title is None:
            title = f"Boxplot: {gene} by {variable}"
        fig = px.box(
            plot_df, x=variable, y=gene, color=variable,
            color_discrete_sequence=px.colors.qualitative.D3,
            color_discrete_map=cmap, range_y=range_y, points=False,
            title=title, template="plotly_white",
        )
        fig.update_xaxes(tickangle=-90, automargin=True)
        fig.update_traces(line_width=1.2, notched=False)
        fig.update_layout(
            xaxis_title=variable, yaxis_title=gene, legend_title=variable,
            width=width, height=height,
        )
    else:
        if adata.isbacked:
            use_adata = adata[:, gene].to_memory()
            adata.file.close()
        else:
            use_adata = adata[:, gene].copy()
        keys = ["min", "q25", "q50", "q75", "max", "mean", "std"]
        plot_data = pd.concat(
            [use_adata.to_df(layer=k)[gene].rename(k) for k in keys], axis=1
        )
        groups = plot_data.sort_values("q50").index.tolist()
        plot_data = plot_data.loc[groups]
        cmap = load_color_palette(palette=palette, adata=adata, groups=variable)
        d3 = px.colors.qualitative.D3
        fig = go.Figure()
        for i, group in enumerate(groups):
            row = plot_data.loc[group]
            color = (cmap or {}).get(group, d3[i % len(d3)])
            fig.add_trace(go.Box(
                x=[group], q1=[row["q25"]], median=[row["q50"]], q3=[row["q75"]],
                lowerfence=[row["min"]], upperfence=[row["max"]],
                boxpoints=False, marker=dict(color=color),
                name=str(group), showlegend=True,
            ))
        if title is None:
            title = f"Boxplot: {gene} by {variable}"
        fig.update_xaxes(tickangle=-90, automargin=True)
        fig.update_layout(
            title=title, xaxis_title=variable, yaxis_title=gene,
            legend_title=variable, template="plotly_white",
            width=width, height=height,
        )
    if show:
        show_fig(fig, filename=f"boxplot.{variable}.{gene}")
        return None
    return fig


# ---------------------------------------------------------------------------
# Multi-gene interactive violin
# ---------------------------------------------------------------------------


def interactive_violin(
    adata,
    genes,
    groupby,
    layer=None,
    use_raw=False,
    palette=None,
    title=None,
    width=1100,
    height=600,
    show=True,
    renderer="notebook",
    box=True,
    points=False,
):
    """Interactive multi-gene violin plot via Plotly.

    Each gene becomes a facet column. Categories along the X axis are
    coloured by ``groupby``.
    """
    import plotly.express as px
    import plotly.io as pio
    import anndata as _ad

    from ._utils import show_fig
    from ..tools import load_color_palette

    if renderer is not None:
        pio.renderers.default = renderer
    if isinstance(genes, str):
        genes = [genes]
    if not isinstance(adata, _ad.AnnData):
        raise TypeError("interactive_violin requires AnnData input.")
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
    cmap = load_color_palette(palette=palette, adata=adata, groups=groupby) if isinstance(
        palette, str
    ) else (palette if isinstance(palette, dict) else None)
    fig = px.violin(
        long, x=groupby, y="value", color=groupby, facet_col="gene",
        facet_col_wrap=min(3, len(genes)), color_discrete_map=cmap,
        color_discrete_sequence=px.colors.qualitative.D3,
        box=box, points="all" if points else False,
        template="plotly_white",
    )
    fig.update_xaxes(tickangle=-45, automargin=True)
    fig.update_yaxes(matches=None, showticklabels=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    if title is None:
        title = f"Expression by {groupby}"
    fig.update_layout(
        title=title, width=width, height=height, legend_title=groupby,
    )
    if show:
        show_fig(fig, filename=f"violin.{groupby}")
        return None
    return fig


# ---------------------------------------------------------------------------
# Stacked bar (composition)
# ---------------------------------------------------------------------------


def interactive_stacked_bar(
    adata,
    groupby,
    split_by,
    normalize=True,
    palette=None,
    order=None,
    split_order=None,
    title=None,
    width=900,
    height=600,
    show=True,
    renderer="notebook",
):
    """Interactive Plotly stacked bar plot of composition."""
    import plotly.express as px
    import plotly.io as pio

    from ._utils import show_fig
    from ..tools import load_color_palette

    if renderer is not None:
        pio.renderers.default = renderer
    obs, ad = resolve_adata_obs(adata)
    splits = categorical_order(obs[split_by], split_order)
    cats = categorical_order(obs[groupby], order)
    ct = pd.crosstab(obs[split_by].astype(str), obs[groupby].astype(str))
    ct = ct.reindex(index=splits, columns=cats, fill_value=0).astype(float)
    if normalize:
        ct = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    long = ct.reset_index().melt(id_vars=split_by, var_name=groupby, value_name="value")
    if isinstance(palette, str):
        cmap = load_color_palette(palette=palette, adata=ad, groups=groupby)
    elif isinstance(palette, dict):
        cmap = palette
    else:
        cmap = None
    fig = px.bar(
        long, x=split_by, y="value", color=groupby,
        category_orders={split_by: splits, groupby: cats},
        color_discrete_map=cmap,
        color_discrete_sequence=px.colors.qualitative.D3,
        template="plotly_white",
    )
    fig.update_xaxes(tickangle=-45, automargin=True)
    fig.update_layout(
        title=title or f"{groupby} composition by {split_by}",
        xaxis_title=split_by,
        yaxis_title="Fraction" if normalize else "Count",
        width=width, height=height, barmode="stack", legend_title=groupby,
    )
    if show:
        show_fig(fig, filename=f"stacked_bar.{split_by}.{groupby}")
        return None
    return fig
