"""Interactive (Plotly) versions of the static cell-statistics plots.

These functions return :class:`plotly.graph_objects.Figure` instances
that can be displayed inline in Jupyter or saved to HTML. They mirror
the static API in :mod:`adataviz.pl.composition` and
:mod:`adataviz.pl.flow`.

All functions lazy-import :mod:`plotly` so that ``import adataviz`` does
not require it.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence, Union

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
    width: Optional[int] = None,
):
    """Interactive polar bar chart (Nightingale rose) of category counts.

    Draws one wedge per category of ``groupby`` whose radius encodes the
    number of cells in that category. Returns an interactive Plotly
    figure that can be displayed inline in Jupyter or written to HTML.

    Parameters
    ----------
    adata : AnnData, DataFrame or path
        Source of the observation metadata. Resolved via
        :func:`resolve_adata_obs`, so an ``AnnData``, an ``obs``-like
        ``DataFrame`` or a path to either is accepted.
    groupby : str
        Column in ``obs`` whose category counts are plotted as the
        angular wedges of the rose.
    order : sequence, optional
        Explicit ordering of the categories around the circle. When
        ``None`` the natural categorical order (or sorted unique values)
        is used.
    palette : dict, str or None, optional
        Colour mapping for the categories. Accepts a ``{category: colour}``
        mapping, a path to an ``.xlsx`` palette sheet, or ``None`` to fall
        back to colours stored in ``adata.uns`` or a generated default.
    title : str, optional
        Figure title. ``None`` leaves the plot untitled.
    save : str, optional
        Output path. A ``.html`` extension writes a standalone HTML file;
        any other extension writes a static image via ``write_image``.
        ``None`` disables saving.
    height : int, default 500
        Figure height in pixels.
    width : int, optional
        Figure width in pixels. ``None`` lets Plotly size it
        automatically.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive polar bar chart figure.
    """
    import plotly.graph_objects as go

    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
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
        width=width,
        margin=dict(l=60, r=60, t=50, b=40),
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
    width: Optional[int] = None,
):
    """Interactive donut chart of ``groupby`` proportions.

    Shows the relative share of each category of ``groupby`` as a pie
    with a central hole. Returns an interactive Plotly figure suitable
    for inline display or HTML export.

    Parameters
    ----------
    adata : AnnData, DataFrame or path
        Source of the observation metadata, resolved via
        :func:`resolve_adata_obs`.
    groupby : str
        Column in ``obs`` whose category proportions are shown as pie
        slices.
    order : sequence, optional
        Explicit ordering of the slices. When ``None`` the natural
        categorical order (or sorted unique values) is used.
    palette : dict, str or None, optional
        Colour mapping for the categories. Accepts a ``{category: colour}``
        mapping, a path to an ``.xlsx`` palette sheet, or ``None`` to fall
        back to colours stored in ``adata.uns`` or a generated default.
    hole : float, default 0.45
        Fraction of the radius cut out at the centre; ``0`` yields a
        full pie and larger values a thinner ring.
    title : str, optional
        Figure title. ``None`` leaves the plot untitled.
    save : str, optional
        Output path. ``.html`` writes standalone HTML; any other
        extension writes a static image. ``None`` disables saving.
    height : int, default 500
        Figure height in pixels.
    width : int, optional
        Figure width in pixels. ``None`` lets Plotly size it
        automatically.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive donut chart figure.
    """
    import plotly.graph_objects as go

    obs, ad = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    counts = obs[groupby].astype(str).value_counts().reindex(cats, fill_value=0)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )
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
    fig.update_layout(
        title=title, height=height, width=width, margin=dict(l=60, r=60, t=40, b=20)
    )
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
    width: Optional[int] = None,
):
    """Interactive pie chart of ``groupby`` proportions (no hole).

    Thin wrapper around :func:`interactive_ring` with ``hole=0.0``, so
    the slices form a solid pie rather than a donut. Returns an
    interactive Plotly figure.

    Parameters
    ----------
    adata : AnnData, DataFrame or path
        Source of the observation metadata, resolved via
        :func:`resolve_adata_obs`.
    groupby : str
        Column in ``obs`` whose category proportions are shown as pie
        slices.
    order : sequence, optional
        Explicit ordering of the slices. When ``None`` the natural
        categorical order (or sorted unique values) is used.
    palette : dict, str or None, optional
        Colour mapping for the categories, forwarded unchanged to
        :func:`interactive_ring`.
    title : str, optional
        Figure title. ``None`` leaves the plot untitled.
    save : str, optional
        Output path. ``.html`` writes standalone HTML; any other
        extension writes a static image. ``None`` disables saving.
    height : int, default 500
        Figure height in pixels.
    width : int, optional
        Figure width in pixels. ``None`` lets Plotly size it
        automatically.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive pie chart figure produced by
        :func:`interactive_ring`.
    """
    return interactive_ring(
        adata,
        groupby=groupby,
        order=order,
        palette=palette,
        hole=0.0,
        title=title,
        save=save,
        height=height,
        width=width,
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
    width: Optional[int] = None,
):
    """Interactive stacked area chart, mirroring :func:`area_plot`.

    Cross-tabulates ``groupby`` against ``split_by`` and stacks one area
    trace per ``groupby`` category across the ``split_by`` values on the
    x-axis. Returns an interactive Plotly figure.

    Parameters
    ----------
    adata : AnnData, DataFrame or path
        Source of the observation metadata, resolved via
        :func:`resolve_adata_obs`.
    groupby : str
        Column in ``obs`` used for the stacked areas (one filled series
        per category).
    split_by : str
        Column in ``obs`` mapped to the x-axis; each value defines a
        column of the underlying cross-tabulation.
    order : sequence, optional
        Explicit ordering of the ``groupby`` categories (stack order).
        ``None`` uses the natural categorical order.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` values along the x-axis.
        ``None`` uses the natural categorical order.
    normalize : bool, default True
        When ``True`` each x-position is normalised to fractions summing
        to 1; when ``False`` raw counts are stacked.
    palette : dict, str or None, optional
        Colour mapping for the ``groupby`` categories. Accepts a
        ``{category: colour}`` mapping, a path to an ``.xlsx`` palette
        sheet, or ``None`` for colours from ``adata.uns`` or a default.
    title : str, optional
        Figure title. ``None`` leaves the plot untitled.
    save : str, optional
        Output path. ``.html`` writes standalone HTML; any other
        extension writes a static image. ``None`` disables saving.
    height : int, default 500
        Figure height in pixels.
    width : int, optional
        Figure width in pixels. ``None`` lets Plotly size it
        automatically.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive stacked area chart figure.
    """
    import plotly.graph_objects as go

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
        width=width,
        xaxis_title=split_by,
        yaxis_title="Fraction" if normalize else "Count",
        margin=dict(l=60, r=20, t=40, b=80),
    )
    fig.update_xaxes(tickangle=-45, automargin=True)
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
    width: Optional[int] = None,
):
    """Interactive line-trend chart, mirroring :func:`trend_plot`.

    Cross-tabulates ``groupby`` against ``split_by`` and draws one
    line-plus-marker trace per ``groupby`` category across the
    ``split_by`` values on the x-axis. Returns an interactive Plotly
    figure.

    Parameters
    ----------
    adata : AnnData, DataFrame or path
        Source of the observation metadata, resolved via
        :func:`resolve_adata_obs`.
    groupby : str
        Column in ``obs`` used for the individual trend lines (one line
        per category).
    split_by : str
        Column in ``obs`` mapped to the x-axis; each value defines a
        column of the underlying cross-tabulation.
    order : sequence, optional
        Explicit ordering of the ``groupby`` categories. ``None`` uses
        the natural categorical order.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` values along the x-axis.
        ``None`` uses the natural categorical order.
    normalize : bool, default True
        When ``True`` each x-position is normalised to fractions summing
        to 1; when ``False`` raw counts are plotted.
    palette : dict, str or None, optional
        Colour mapping for the ``groupby`` categories. Accepts a
        ``{category: colour}`` mapping, a path to an ``.xlsx`` palette
        sheet, or ``None`` for colours from ``adata.uns`` or a default.
    title : str, optional
        Figure title. ``None`` leaves the plot untitled.
    save : str, optional
        Output path. ``.html`` writes standalone HTML; any other
        extension writes a static image. ``None`` disables saving.
    height : int, default 500
        Figure height in pixels.
    width : int, optional
        Figure width in pixels. ``None`` lets Plotly size it
        automatically.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive line-trend chart figure.
    """
    import plotly.graph_objects as go

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
        width=width,
        xaxis_title=split_by,
        yaxis_title="Fraction" if normalize else "Count",
        margin=dict(l=60, r=20, t=40, b=80),
    )
    fig.update_xaxes(tickangle=-45, automargin=True)
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
    width: Optional[int] = None,
):
    """Interactive Sankey diagram mirroring :func:`sankey_plot`.

    Cross-tabulates the ``left`` and ``right`` columns and draws a flow
    diagram whose link widths encode the number of cells transitioning
    between each pair of categories; empty nodes are dropped and each
    link can be hovered for its underlying count. Returns an interactive
    Plotly figure.

    Parameters
    ----------
    adata : AnnData, DataFrame or path
        Source of the observation metadata, resolved via
        :func:`resolve_adata_obs`.
    left : str
        Column in ``obs`` providing the left-hand (source) nodes.
    right : str
        Column in ``obs`` providing the right-hand (target) nodes.
    left_order : sequence, optional
        Explicit ordering of the left-hand nodes. ``None`` uses the
        natural categorical order.
    right_order : sequence, optional
        Explicit ordering of the right-hand nodes. ``None`` uses the
        natural categorical order.
    palette : dict, str or None, optional
        Colour mapping for the left-hand nodes (link colours are derived
        from these as semi-transparent fills); right-hand nodes are
        coloured from ``adata.uns`` or a default. Accepts a
        ``{category: colour}`` mapping, an ``.xlsx`` palette path, or
        ``None``.
    title : str, optional
        Figure title. ``None`` leaves the plot untitled.
    save : str, optional
        Output path. ``.html`` writes standalone HTML; any other
        extension writes a static image. ``None`` disables saving.
    height : int, default 500
        Figure height in pixels.
    width : int, optional
        Figure width in pixels. ``None`` lets Plotly size it
        automatically.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive Sankey diagram figure.
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

    left_colors = resolve_palette(
        palette, left_cats, sheet_name=left, adata=ad, groupby=left
    )
    right_colors = resolve_palette(
        None, right_cats, sheet_name=right, adata=ad, groupby=right
    )
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
        title=title, height=height, width=width, margin=dict(l=60, r=60, t=40, b=20)
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
    """Interactive scatter on a 2-D embedding (UMAP, t-SNE, etc.).

    Colours points either by a categorical ``obs`` column (``variable``)
    or by a continuous gene's expression (``gene``). The figure is built
    with :func:`plotly.express.scatter` using a WebGL renderer so it can
    handle large numbers of cells interactively.

    Parameters
    ----------
    adata : AnnData or path, optional
        Source data. Required when ``gene`` is supplied (expression is
        read from it) and used for the embedding coordinates; may be
        ``None`` if ``obs`` already contains the embedding columns.
    obs : DataFrame or path, optional
        Pre-computed observation table. When given it is intersected with
        ``adata`` (if any); otherwise ``adata.obs`` is used.
    variable : str, optional
        Categorical ``obs`` column to colour by. Mutually intended with
        ``gene`` (``gene`` takes precedence when both are set).
    gene : str, optional
        Gene name to colour by continuous expression; requires
        ``adata``.
    annotation : str, optional
        Extra ``obs`` column added to the hover tooltip only.
    basis : str, default "umap"
        Embedding key (without the ``X_`` prefix) read from
        ``adata.obsm[f"X_{basis}"]``.
    vmin : str, default "p1"
        Lower percentile (``"p<N>"``) used as the low end of the
        continuous colour range for ``gene`` colouring.
    vmax : str, default "p99"
        Upper percentile (``"p<N>"``) used as the high end of the
        continuous colour range for ``gene`` colouring.
    cmap : str, default "jet"
        Continuous colour scale name passed to Plotly for gene
        expression.
    title : str, optional
        Figure title. ``None`` auto-generates one from ``basis`` and the
        coloured field.
    width : int, default 900
        Figure width in pixels (also feeds automatic marker sizing).
    height : int, default 750
        Figure height in pixels (also feeds automatic marker sizing).
    palette : dict, str or None, optional
        Categorical colour mapping. A ``{category: colour}`` dict or a
        path to an ``.xlsx`` palette sheet (resolved via
        :func:`load_color_palette`); unused categories are dropped.
        ``None`` falls back to a qualitative default sequence.
    size : int, optional
        Marker diameter in pixels. ``None`` derives a size automatically
        from the point count and figure area.
    show : bool, default True
        When ``True`` the figure is displayed via ``show_fig`` and
        ``None`` is returned; when ``False`` the figure object is
        returned instead.
    downsample : int, optional
        If set and the data has more rows than this, randomly subsample
        to this many points before plotting.
    target_fill : float, default 0.05
        Target fraction of the plot area covered by markers, used by the
        automatic marker-size heuristic when ``size`` is ``None``.
    normalize_per_cell : bool, default False
        When colouring by ``gene``, normalise expression per cell (via
        :func:`normalize_mc_by_cell`) before plotting.
    clip_norm_value : float, default 10
        Clipping ceiling applied during per-cell normalisation.
    renderer : str, default "notebook"
        Plotly IO renderer to activate (e.g. ``"notebook"``,
        ``"browser"``). ``None`` leaves the current default unchanged.

    Returns
    -------
    plotly.graph_objects.Figure or None
        The interactive embedding scatter figure, or ``None`` when
        ``show=True`` (the figure is displayed instead of returned).
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
                use_adata=use_adata,
                normalize_per_cell=normalize_per_cell,
                clip_norm_value=clip_norm_value,
                hypo_score=False,
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
        obs,
        x=f"{basis}_0",
        y=f"{basis}_1",
        color=use_col,
        category_orders=category_orders,
        hover_data=hover_data,
        range_color=range_color,
        color_discrete_sequence=px.colors.qualitative.D3,
        color_discrete_map=color_discrete_map,
        color_continuous_scale=cmap,
        template="plotly_white",
        render_mode="webgl",
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
        marker_diam_area = 2 * np.sqrt(
            (width * height * target_fill) / (np.pi * n_points)
        )
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
        autosize=True,
        width=width,
        height=height,
        legend_title=use_col,
        legend=dict(font_size=12, itemsizing="constant", itemwidth=30, borderwidth=0.1),
    )
    if show:
        show_fig(fig, filename=f"{basis}.{use_col}")
        return None
    return fig


def interactive_categorical(adata, groupby, basis="umap", **kwargs):
    """Interactive embedding scatter coloured by a categorical column.

    Thin convenience wrapper around :func:`interactive_embedding` that
    forwards ``groupby`` as the categorical ``variable`` argument.

    Parameters
    ----------
    adata : AnnData or path
        Source data providing the embedding and ``obs``.
    groupby : str
        Categorical ``obs`` column to colour by (passed as ``variable``).
    basis : str, default "umap"
        Embedding key (without the ``X_`` prefix).
    **kwargs
        Additional keyword arguments forwarded verbatim to
        :func:`interactive_embedding` (e.g. ``palette``, ``size``,
        ``show``, ``width``, ``height``, ``renderer``).

    Returns
    -------
    plotly.graph_objects.Figure or None
        The figure from :func:`interactive_embedding`, or ``None`` when
        it is asked to display the plot (``show=True``).
    """
    return interactive_embedding(adata=adata, variable=groupby, basis=basis, **kwargs)


def interactive_gene(adata, gene, basis="umap", **kwargs):
    """Interactive embedding scatter coloured by a gene's expression.

    Thin convenience wrapper around :func:`interactive_embedding` that
    forwards ``gene`` for continuous expression colouring.

    Parameters
    ----------
    adata : AnnData or path
        Source data; must contain ``gene`` in ``var_names`` and the
        requested embedding.
    gene : str
        Gene name to colour by continuous expression.
    basis : str, default "umap"
        Embedding key (without the ``X_`` prefix).
    **kwargs
        Additional keyword arguments forwarded verbatim to
        :func:`interactive_embedding` (e.g. ``cmap``, ``vmin``, ``vmax``,
        ``normalize_per_cell``, ``show``, ``renderer``).

    Returns
    -------
    plotly.graph_objects.Figure or None
        The figure from :func:`interactive_embedding`, or ``None`` when
        it is asked to display the plot (``show=True``).
    """
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
    query=None,
):
    """Interactive dot heatmap of gene expression across groups.

    For each ``(group, gene)`` pair the dot size encodes the fraction of
    expressing cells and the dot colour encodes the mean expression.
    Expression statistics are prepared by :func:`get_genes_mean_frac`
    (which also handles normalisation and the ``query`` subset), and the
    result is an interactive Plotly figure.

    Parameters
    ----------
    adata : AnnData or path, optional
        Source data passed to :func:`get_genes_mean_frac`. May be
        ``None`` if ``obs`` supplies the needed metadata.
    obs : DataFrame or path, optional
        Alternative observation metadata forwarded to
        :func:`get_genes_mean_frac`.
    genes : sequence of str, optional
        Genes to include. ``None`` defers gene selection to
        :func:`get_genes_mean_frac`.
    groupby : str, default "Subclass"
        ``obs`` column defining the groups placed along the x-axis.
    modality : str, default "RNA"
        Data modality forwarded to :func:`get_genes_mean_frac` for
        selecting the expression matrix.
    title : str, optional
        Figure title. ``None`` defaults to ``groupby``.
    use_raw : bool, default False
        Whether to use ``adata.raw`` when computing expression, forwarded
        to :func:`get_genes_mean_frac`.
    expression_cutoff : str, default "p5"
        Percentile-style cutoff (``"p<N>"``) used by
        :func:`get_genes_mean_frac` to decide which cells count as
        expressing.
    normalize_per_cell : bool, default True
        Whether to normalise expression per cell before summarising,
        forwarded to :func:`get_genes_mean_frac`.
    clip_norm_value : float, default 10
        Clipping ceiling applied during per-cell normalisation.
    width : int, default 900
        Minimum figure width in pixels; the actual width grows with the
        number of groups so dots and labels stay legible.
    height : int, default 700
        Minimum figure height in pixels; the actual height grows with the
        number of genes.
    gene_order : sequence of str, optional
        Explicit ordering of genes on the y-axis (only genes present in
        the data are kept). ``None`` uses the order returned by
        :func:`get_genes_mean_frac`.
    colorscale : str, default "greens"
        Plotly colour scale used to encode mean expression.
    vmin : str, default "p1"
        Lower percentile (``"p<N>"``) of mean expression used as the low
        end of the colour range.
    vmax : str, default "p99"
        Upper percentile (``"p<N>"``) of mean expression used as the high
        end of the colour range.
    show : bool, default True
        When ``True`` the figure is displayed via ``show_fig`` and
        ``None`` is returned; when ``False`` the figure object is
        returned.
    reversescale : bool, default False
        Reverse the colour scale direction.
    size_min : int, default 5
        Marker size (pixels) for a fraction of 0.
    size_max : int, default 30
        Marker size (pixels) for a fraction of 1; intermediate fractions
        are linearly interpolated between ``size_min`` and ``size_max``.
    renderer : str, default "notebook"
        Plotly IO renderer to activate. ``None`` leaves the current
        default unchanged.
    query : str, optional
        Pandas-style query string passed through to
        :func:`get_genes_mean_frac` to subset ``adata.obs`` before
        materialising the expression matrix. Example:
        ``"Subclass == 'L2/3'"``.

    Returns
    -------
    plotly.graph_objects.Figure or None
        The interactive dot heatmap figure, or ``None`` when
        ``show=True`` (the figure is displayed instead of returned).
    """
    import plotly.io as pio
    import plotly.graph_objects as go

    from ._utils import show_fig
    from .expression import get_genes_mean_frac

    if renderer is not None:
        pio.renderers.default = renderer

    plot_data = get_genes_mean_frac(
        adata,
        obs=obs,
        groupby=groupby,
        modality=modality,
        use_raw=use_raw,
        expression_cutoff=expression_cutoff,
        genes=genes,
        normalize_per_cell=normalize_per_cell,
        clip_norm_value=clip_norm_value,
        hypo_score=False,
        query=query,
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
            plot_data[groupby].tolist(),
            plot_data["Gene"].tolist(),
            mean_vals,
            frac_vals,
        )
    ]
    vmin_q = float(int(vmin.replace("p", "")) / 100)
    vmax_q = float(int(vmax.replace("p", "")) / 100)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=plot_data[groupby].tolist(),
            y=plot_data["Gene"].tolist(),
            mode="markers",
            marker=dict(
                size=sizes,
                color=mean_vals,
                colorscale=colorscale,
                showscale=True,
                colorbar=dict(title="Mean"),
                reversescale=reversescale,
                sizemode="area",
                opacity=0.9,
                cmin=plot_data["Mean"].quantile(vmin_q),
                cmax=plot_data["Mean"].quantile(vmax_q),
            ),
            text=hover_text,
            hoverinfo="text",
        )
    )
    fig.update_xaxes(
        type="category", categoryorder="array", categoryarray=x_labels,
        tickangle=-45, automargin=True,
    )
    fig.update_yaxes(
        type="category", categoryorder="array",
        categoryarray=list(reversed(y_labels)), automargin=True,
    )
    if title is None:
        title = groupby
    # Grow with the number of groups (x) and genes (y) so dots/labels don't
    # get crushed; reserve margins for angled labels and the colorbar.
    auto_w = max(width, 250 + 28 * len(x_labels))
    auto_h = max(height, 150 + 24 * len(y_labels))
    fig.update_layout(
        title=title or "",
        xaxis_title=groupby,
        yaxis_title="Gene",
        width=auto_w,
        height=auto_h,
        plot_bgcolor="white",
        margin=dict(l=80, r=90, t=60, b=120),
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
    return all(
        k in adata.layers for k in ["min", "q25", "q50", "q75", "max", "mean", "std"]
    )


def _get_boxplot_data(adata, variable, gene, obs=None):
    import os
    import anndata

    assert isinstance(adata, anndata.AnnData)
    if isinstance(obs, str):
        obs_path = os.path.abspath(os.path.expanduser(obs))
        sep = "\t" if obs_path.endswith((".tsv", ".txt")) else ","
        data = pd.read_csv(obs_path, index_col=0, sep=sep)
    else:
        assert isinstance(obs, pd.DataFrame)
        data = obs.copy()
    overlap = data.index.intersection(adata.obs_names)
    data = data.loc[overlap]
    # Subset cells + single gene before materialising, so backed AnnData
    # only loads `len(overlap)` rows. h5py disallows fancy indexing on
    # both axes at once, so subset rows on the backed object first and
    # take the column after materialisation.
    if adata.isbacked:
        use_adata = adata[overlap].to_memory()[:, gene].copy()
    else:
        use_adata = adata[overlap, gene].copy()
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
    query=None,
):
    """Interactive Plotly boxplot of a single gene grouped by a categorical variable.

    Builds an interactive box plot for one ``gene`` split across the
    categories of ``variable``. The function auto-detects whether
    ``adata`` carries pre-computed statistic layers (``min``/``q25``/
    ``q50``/``q75``/``max``/``mean``/``std``, i.e. pseudobulk); if so it
    draws pre-summarised boxes, otherwise it computes the box from the
    raw per-cell values. For multi-gene support use
    :func:`interactive_violin`.

    Parameters
    ----------
    adata : anndata.AnnData or str
        Annotated data matrix, or a path to an ``.h5ad`` file. A string
        is opened in backed (``"r"``) mode so that only the required
        rows/columns are materialised.
    variable : str
        Column in ``adata.obs`` (or ``obs``) whose categories define the
        X-axis groups and the box colours.
    gene : str
        Gene (``adata.var_names`` entry) plotted on the Y axis.
    obs : pandas.DataFrame, optional
        Pre-built observation table to use instead of ``adata.obs``. When
        ``None`` (default) the observations are taken from ``adata``
        (subset by ``query`` if given).
    palette : str or dict, optional
        Colour specification passed to
        :func:`adataviz.tools.load_color_palette`. May be a palette name
        resolved against ``adata`` for ``variable`` or a ``{category:
        colour}`` mapping. When ``None`` a default qualitative (D3)
        sequence is used.
    title : str, optional
        Figure title. When ``None`` a title of the form
        ``"Boxplot: {gene} by {variable}"`` is generated.
    width : int, default 1100
        Figure width in pixels.
    height : int, default 700
        Figure height in pixels.
    show : bool, default True
        If ``True`` the figure is rendered via the internal ``show_fig``
        helper and ``None`` is returned; if ``False`` the
        :class:`plotly.graph_objects.Figure` is returned instead.
    renderer : str, default "notebook"
        Plotly renderer to set as the default (``plotly.io.renderers``).
        Pass ``None`` to leave the current renderer unchanged.
    query : str, optional
        Pandas-style query string applied to ``adata.obs`` (and to
        ``obs`` if provided) to subset cells before computing the
        boxplot. For backed AnnData this avoids materialising rows
        that will be discarded. Example: ``"Region == 'CTX'"``.

    Returns
    -------
    plotly.graph_objects.Figure or None
        The interactive box-plot figure, or ``None`` when ``show`` is
        ``True`` (the figure is displayed instead of returned).
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
    # Resolve query into a kept_cells index, but defer the actual row
    # selection until the per-gene to_memory step below so that backed
    # AnnData only loads `len(kept_cells)` rows of a single column.
    kept_cells = None
    if query is not None:
        kept_cells = adata.obs.query(query).index
        if obs is not None and isinstance(obs, pd.DataFrame):
            obs = obs.query(query)
    if obs is None:
        obs = (
            adata.obs.loc[kept_cells].copy()
            if kept_cells is not None
            else adata.obs.copy()
        )
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
            plot_df,
            x=variable,
            y=gene,
            color=variable,
            color_discrete_sequence=px.colors.qualitative.D3,
            color_discrete_map=cmap,
            range_y=range_y,
            points=False,
            title=title,
            template="plotly_white",
        )
        fig.update_xaxes(tickangle=-90, automargin=True)
        fig.update_traces(line_width=1.2, notched=False)
        fig.update_layout(
            xaxis_title=variable,
            yaxis_title=gene,
            legend_title=variable,
            width=width,
            height=height,
        )
    else:
        row_sel = kept_cells if kept_cells is not None else slice(None)
        if adata.isbacked:
            # h5py disallows fancy indexing on both axes simultaneously.
            use_adata = adata[row_sel].to_memory()[:, gene].copy()
            adata.file.close()
        else:
            use_adata = adata[row_sel, gene].copy()
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
            fig.add_trace(
                go.Box(
                    x=[group],
                    q1=[row["q25"]],
                    median=[row["q50"]],
                    q3=[row["q75"]],
                    lowerfence=[row["min"]],
                    upperfence=[row["max"]],
                    boxpoints=False,
                    marker=dict(color=color),
                    name=str(group),
                    showlegend=True,
                )
            )
        if title is None:
            title = f"Boxplot: {gene} by {variable}"
        fig.update_xaxes(tickangle=-90, automargin=True)
        fig.update_layout(
            title=title,
            xaxis_title=variable,
            yaxis_title=gene,
            legend_title=variable,
            template="plotly_white",
            width=width,
            height=height,
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

    Each gene becomes a facet column and the categories of ``groupby``
    are laid out along the X axis and coloured accordingly. Expression
    values are melted into long form so a single faceted figure spans
    all requested genes.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix. A concrete ``AnnData`` is required (a
        :class:`TypeError` is raised otherwise).
    genes : str or sequence of str
        Gene or genes to plot. A single string is wrapped into a
        one-element list; each gene becomes a facet column.
    groupby : str
        Column in ``adata.obs`` defining the X-axis categories and the
        violin colours.
    layer : str, optional
        Name of an ``adata.layers`` entry to read expression from. When
        ``None`` (default) the main ``X`` matrix (or ``raw.X`` if
        ``use_raw``) is used.
    use_raw : bool, default False
        If ``True`` and ``adata.raw`` is set, read gene values and
        ``var_names`` from ``adata.raw`` instead of ``adata``.
    palette : str or dict, optional
        Colour specification. A string is resolved via
        :func:`adataviz.tools.load_color_palette` for ``groupby``; a
        dict is used directly as a ``{category: colour}`` map. When
        ``None`` a default qualitative (D3) sequence is used.
    title : str, optional
        Figure title. When ``None`` a title of the form
        ``"Expression by {groupby}"`` is generated.
    width : int, default 1100
        Figure width in pixels.
    height : int, default 600
        Figure height in pixels.
    show : bool, default True
        If ``True`` the figure is displayed via ``show_fig`` and ``None``
        is returned; if ``False`` the figure object is returned.
    renderer : str, default "notebook"
        Plotly renderer to set as the default. Pass ``None`` to leave the
        current renderer unchanged.
    box : bool, default True
        Whether to overlay a box plot inside each violin.
    points : bool, default False
        If ``True`` all underlying data points are shown (``points="all"``);
        if ``False`` no points are drawn.

    Returns
    -------
    plotly.graph_objects.Figure or None
        The interactive faceted violin figure, or ``None`` when ``show``
        is ``True``.
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
    cmap = (
        load_color_palette(palette=palette, adata=adata, groups=groupby)
        if isinstance(palette, str)
        else (palette if isinstance(palette, dict) else None)
    )
    fig = px.violin(
        long,
        x=groupby,
        y="value",
        color=groupby,
        facet_col="gene",
        facet_col_wrap=min(3, len(genes)),
        color_discrete_map=cmap,
        color_discrete_sequence=px.colors.qualitative.D3,
        box=box,
        points="all" if points else False,
        template="plotly_white",
    )
    fig.update_xaxes(tickangle=-45, automargin=True)
    fig.update_yaxes(matches=None, showticklabels=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    if title is None:
        title = f"Expression by {groupby}"
    fig.update_layout(
        title=title,
        width=width,
        height=height,
        legend_title=groupby,
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
    """Interactive Plotly stacked bar plot of categorical composition.

    Cross-tabulates ``groupby`` against ``split_by`` and draws one
    stacked bar per ``split_by`` category, with each stack segment
    representing a ``groupby`` category. Optionally normalises each bar
    to fractions summing to one.

    Parameters
    ----------
    adata : anndata.AnnData or str
        Annotated data matrix or path; the observation table is resolved
        via ``resolve_adata_obs``.
    groupby : str
        Column in ``adata.obs`` whose categories become the coloured
        stack segments within each bar.
    split_by : str
        Column in ``adata.obs`` whose categories define the individual
        bars along the X axis.
    normalize : bool, default True
        If ``True`` each bar is normalised so segment heights are
        fractions summing to one (Y axis labelled "Fraction"); if
        ``False`` raw counts are shown (Y axis labelled "Count").
    palette : str or dict, optional
        Colour specification for ``groupby``. A string is resolved via
        :func:`adataviz.tools.load_color_palette`; a dict is used as a
        ``{category: colour}`` map. When ``None`` a default qualitative
        (D3) sequence is used.
    order : sequence, optional
        Explicit ordering of the ``groupby`` categories (stack order and
        legend). When ``None`` the natural categorical order is used.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` categories (bar order along
        the X axis). When ``None`` the natural categorical order is used.
    title : str, optional
        Figure title. When ``None`` a title of the form
        ``"{groupby} composition by {split_by}"`` is generated.
    width : int, default 900
        Figure width in pixels.
    height : int, default 600
        Figure height in pixels.
    show : bool, default True
        If ``True`` the figure is displayed via ``show_fig`` and ``None``
        is returned; if ``False`` the figure object is returned.
    renderer : str, default "notebook"
        Plotly renderer to set as the default. Pass ``None`` to leave the
        current renderer unchanged.

    Returns
    -------
    plotly.graph_objects.Figure or None
        The interactive stacked bar figure, or ``None`` when ``show`` is
        ``True``.
    """
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
        long,
        x=split_by,
        y="value",
        color=groupby,
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
        width=width,
        height=height,
        barmode="stack",
        legend_title=groupby,
    )
    if show:
        show_fig(fig, filename=f"stacked_bar.{split_by}.{groupby}")
        return None
    return fig


# ---------------------------------------------------------------------------
# Dot plot (count × column-fraction)
# ---------------------------------------------------------------------------


def interactive_dot_plot(
    adata,
    groupby: str,
    split_by: str,
    order: Optional[Sequence] = None,
    split_order: Optional[Sequence] = None,
    cmap: str = "Viridis",
    size_max: int = 30,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
    width: int = 700,
):
    """Interactive Plotly dot plot where marker size = count and colour = column-fraction.

    Cross-tabulates ``groupby`` (rows) against ``split_by`` (columns) and
    draws a grid of dots: dot area encodes the raw cell count while dot
    colour encodes the fraction of that count within its ``split_by``
    column. Mirrors :func:`adataviz.pl.dot_plot`.

    Parameters
    ----------
    adata : anndata.AnnData or str
        Annotated data matrix or path; the observation table is resolved
        via ``resolve_adata_obs``.
    groupby : str
        Column in ``adata.obs`` whose categories form the Y-axis rows.
    split_by : str
        Column in ``adata.obs`` whose categories form the X-axis columns
        and define the denominator for the colour fractions.
    order : sequence, optional
        Explicit ordering of the ``groupby`` categories (Y axis). When
        ``None`` the natural categorical order is used.
    split_order : sequence, optional
        Explicit ordering of the ``split_by`` categories (X axis). When
        ``None`` the natural categorical order is used.
    cmap : str, default "Viridis"
        Plotly continuous colourscale name used to map column fractions
        (0-1) onto dot colours.
    size_max : int, default 30
        Maximum dot size in pixels; scales the ``sizeref`` used to map
        counts to marker area.
    title : str, optional
        Figure title. When ``None`` no title is shown.
    save : str, optional
        Path to save the figure to via the internal ``_maybe_save``
        helper. When ``None`` the figure is not written to disk.
    height : int, default 500
        Minimum figure height in pixels; the effective height grows with
        the number of ``groupby`` categories.
    width : int, default 700
        Minimum figure width in pixels; the effective width grows with
        the number of ``split_by`` categories.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive dot-plot figure.
    """
    import plotly.graph_objects as go

    obs, _ = resolve_adata_obs(adata)
    cats = categorical_order(obs[groupby], order)
    splits = categorical_order(obs[split_by], split_order)
    ct = pd.crosstab(obs[groupby].astype(str), obs[split_by].astype(str))
    ct = ct.reindex(index=cats, columns=splits, fill_value=0).astype(float)
    col_sums = ct.sum(axis=0).replace(0, np.nan)
    frac = ct.div(col_sums, axis=1).fillna(0)
    long = ct.stack().rename("count").reset_index()
    long.columns = [groupby, split_by, "count"]
    long["fraction"] = frac.stack().values

    fig = go.Figure(
        go.Scatter(
            x=long[split_by],
            y=long[groupby],
            mode="markers",
            marker=dict(
                size=long["count"],
                sizemode="area",
                sizeref=2.0 * max(long["count"].max(), 1) / (size_max**2),
                sizemin=2,
                color=long["fraction"],
                colorscale=cmap,
                cmin=0,
                cmax=1,
                line=dict(width=0.5, color="black"),
                colorbar=dict(title=f"Fraction within {split_by}"),
            ),
            customdata=np.stack([long["count"], long["fraction"]], axis=-1),
            hovertemplate=(
                f"{groupby}: %{{y}}<br>{split_by}: %{{x}}<br>"
                "count: %{customdata[0]}<br>fraction: %{customdata[1]:.3f}<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title=split_by,
        yaxis_title=groupby,
        yaxis=dict(autorange="reversed"),
        height=max(height, 150 + 22 * len(cats)),
        width=max(width, 250 + 30 * len(splits)),
        template="plotly_white",
        margin=dict(l=120, r=90, t=50, b=100),
    )
    fig.update_xaxes(tickangle=-45, automargin=True)
    fig.update_yaxes(automargin=True)
    _maybe_save(fig, save)
    return fig


# ---------------------------------------------------------------------------
# Chord (interactive proxy via Sankey)
# ---------------------------------------------------------------------------


def interactive_chord(
    adata,
    left: str,
    right: str,
    left_order: Optional[Sequence] = None,
    right_order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 600,
    width: int = 800,
):
    """Interactive co-occurrence visualisation between two categorical columns.

    Plotly does not provide a native chord renderer, so this function
    draws a Sankey with both columns as separate node sets - the same
    information a chord plot encodes circularly. It is a thin wrapper that
    forwards all arguments to :func:`interactive_sankey`.

    Parameters
    ----------
    adata : anndata.AnnData or str
        Annotated data matrix or path; forwarded to
        :func:`interactive_sankey`.
    left : str
        Column in ``adata.obs`` used as the left-hand node set (source).
    right : str
        Column in ``adata.obs`` used as the right-hand node set (target).
    left_order : sequence, optional
        Explicit ordering of the ``left`` categories. When ``None`` the
        natural categorical order is used.
    right_order : sequence, optional
        Explicit ordering of the ``right`` categories. When ``None`` the
        natural categorical order is used.
    palette : dict, str or None, optional
        Colour specification for the nodes, forwarded to
        :func:`interactive_sankey`. May be a ``{category: colour}``
        mapping, a palette name, or ``None`` for defaults.
    title : str, optional
        Figure title. When ``None`` the wrapped function's default is
        used.
    save : str, optional
        Path to save the figure to. When ``None`` the figure is not
        written to disk.
    height : int, default 600
        Figure height in pixels.
    width : int, default 800
        Figure width in pixels.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive Sankey figure produced by
        :func:`interactive_sankey`.
    """
    return interactive_sankey(
        adata,
        left=left,
        right=right,
        left_order=left_order,
        right_order=right_order,
        palette=palette,
        title=title,
        save=save,
        height=height,
        width=width,
    )


# ---------------------------------------------------------------------------
# Venn / UpSet (interactive bar-style overlap viewer)
# ---------------------------------------------------------------------------


def interactive_upset(
    adata,
    groupby: str,
    set_by: str,
    order: Optional[Sequence] = None,
    palette: Union[None, Mapping[str, str], str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 600,
    width: int = 900,
    top_n: int = 25,
    min_intersection_size: int = 1,
):
    """Interactive UpSet plot of set overlaps across a categorical variable.

    For each category of ``groupby`` a set is built from the unique
    ``set_by`` values it contains, and the sizes of every observed
    combination of sets are visualised. Renders three coordinated panels
    via :mod:`plotly.subplots`:

    - **Top**: vertical bars with the size of each intersection.
    - **Bottom-right**: a dot/line matrix indicating which sets each
      intersection belongs to (filled dots = "in", grey dots = "out",
      vertical line connects all "in" dots in a column).
    - **Bottom-left**: horizontal bars with the total size of every set,
      coloured by ``palette``.

    Parameters
    ----------
    adata : anndata.AnnData or str
        Annotated data matrix or path; the observation table is resolved
        via ``resolve_adata_obs``.
    groupby : str
        Column in ``adata.obs`` whose categories define the individual
        sets.
    set_by : str
        Column in ``adata.obs`` whose values are the items distributed
        into the sets (the universe of elements whose overlaps are
        counted).
    order : sequence, optional
        Explicit ordering of the ``groupby`` categories (set order). When
        ``None`` the natural categorical order is used.
    palette : dict, str or None, optional
        Colour specification for the sets, resolved via
        ``resolve_palette``. May be a ``{category: colour}`` mapping, a
        palette name, or ``None`` for defaults.
    title : str, optional
        Figure title. When ``None`` a title of the form
        ``"{set_by} overlap across {groupby}"`` is generated.
    save : str, optional
        Path to save the figure to via ``_maybe_save``. When ``None`` the
        figure is not written to disk.
    height : int, default 600
        Figure height in pixels.
    width : int, default 900
        Figure width in pixels.
    top_n : int, default 25
        Keep the ``top_n`` largest intersections.
    min_intersection_size : int, default 1
        Drop combinations seen fewer times than this threshold.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive UpSet figure.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    obs, ad = resolve_adata_obs(adata)
    groups = list(categorical_order(obs[groupby], order))
    sets = {
        g: set(obs.loc[obs[groupby].astype(str) == g, set_by].dropna().astype(str))
        for g in groups
    }
    set_sizes = {g: len(sets[g]) for g in groups}
    universe = set().union(*sets.values())
    if not universe:
        raise ValueError("interactive_upset: empty universe of items.")

    # Compute intersection counts.
    combo_counts: dict = {}
    for item in universe:
        combo = tuple(g for g in groups if item in sets[g])
        if not combo:
            continue
        combo_counts[combo] = combo_counts.get(combo, 0) + 1
    rows = [
        (combo, n) for combo, n in combo_counts.items() if n >= min_intersection_size
    ]
    rows.sort(key=lambda r: r[1], reverse=True)
    rows = rows[:top_n]
    if not rows:
        raise ValueError(
            "interactive_upset: no combinations passed `min_intersection_size`."
        )

    n_int = len(rows)
    n_set = len(groups)
    int_sizes = [n for _, n in rows]
    int_labels = [" & ".join(combo) for combo, _ in rows]
    colors = resolve_palette(
        palette, groups, sheet_name=groupby, adata=ad, groupby=groupby
    )
    set_colors = [colors[g] for g in groups]

    # Layout: 2 rows x 2 cols.
    # row1: empty | top bars
    # row2: side bars | dot matrix
    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.22, 0.78],
        row_heights=[0.55, 0.45],
        horizontal_spacing=0.04,
        vertical_spacing=0.06,
        specs=[
            [None, {"type": "bar"}],
            [{"type": "bar"}, {"type": "scatter"}],
        ],
        shared_xaxes=False,
        shared_yaxes=False,
    )

    # --- Top: intersection-size bars --------------------------------------
    fig.add_trace(
        go.Bar(
            x=list(range(n_int)),
            y=int_sizes,
            marker=dict(color="#444"),
            text=int_sizes,
            textposition="outside",
            hovertemplate="<b>%{customdata}</b><br>size: %{y}<extra></extra>",
            customdata=int_labels,
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    # --- Bottom-left: per-set total size bars ------------------------------
    fig.add_trace(
        go.Bar(
            y=list(range(n_set)),
            x=[set_sizes[g] for g in groups],
            orientation="h",
            marker=dict(color=set_colors),
            text=[set_sizes[g] for g in groups],
            textposition="outside",
            hovertemplate="<b>%{customdata}</b><br>set size: %{x}<extra></extra>",
            customdata=groups,
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    # --- Bottom-right: dot matrix ------------------------------------------
    on_x, on_y = [], []
    off_x, off_y = [], []
    seg_x, seg_y = [], []
    for ci, (combo, _) in enumerate(rows):
        in_idx = [groups.index(g) for g in combo]
        for gi in range(n_set):
            if gi in in_idx:
                on_x.append(ci)
                on_y.append(gi)
            else:
                off_x.append(ci)
                off_y.append(gi)
        if len(in_idx) >= 2:
            seg_x += [ci, ci, None]
            seg_y += [min(in_idx), max(in_idx), None]
    fig.add_trace(
        go.Scatter(
            x=off_x,
            y=off_y,
            mode="markers",
            marker=dict(size=10, color="#dddddd", line=dict(color="white", width=1)),
            hoverinfo="skip",
            showlegend=False,
        ),
        row=2,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=seg_x,
            y=seg_y,
            mode="lines",
            line=dict(color="#222", width=2),
            hoverinfo="skip",
            showlegend=False,
        ),
        row=2,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=on_x,
            y=on_y,
            mode="markers",
            marker=dict(
                size=12,
                color=[set_colors[g] for g in on_y],
                line=dict(color="black", width=0.5),
            ),
            hovertemplate="set: %{text}<extra></extra>",
            text=[groups[g] for g in on_y],
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    fig.update_xaxes(showticklabels=False, range=[-0.5, n_int - 0.5], row=1, col=2)
    # Headroom so the "outside" bar count labels aren't clipped by the title.
    fig.update_yaxes(
        title_text="Intersection size",
        range=[0, max(int_sizes) * 1.15] if int_sizes else None,
        row=1,
        col=2,
    )
    fig.update_xaxes(title_text="Set size", autorange="reversed", row=2, col=1)
    fig.update_yaxes(
        tickmode="array",
        tickvals=list(range(n_set)),
        ticktext=groups,
        row=2,
        col=1,
    )
    fig.update_xaxes(
        showticklabels=False,
        range=[-0.5, n_int - 0.5],
        row=2,
        col=2,
    )
    fig.update_yaxes(
        tickmode="array",
        tickvals=list(range(n_set)),
        ticktext=["" for _ in groups],
        range=[-0.5, n_set - 0.5],
        row=2,
        col=2,
    )
    fig.update_layout(
        title=title or f"{set_by} overlap across {groupby}",
        template="plotly_white",
        height=height,
        width=width,
        margin=dict(l=60, r=40, t=60, b=60),
        bargap=0.25,
    )
    _maybe_save(fig, save)
    return fig


# ---------------------------------------------------------------------------
# Complex heatmap (Plotly heatmap of mean expression)
# ---------------------------------------------------------------------------


def interactive_complex_heatmap(
    adata,
    genes: Sequence[str],
    groupby: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    z_score: Optional[str] = "row",
    cmap: str = "RdBu_r",
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
    width: int = 800,
):
    """Interactive Plotly heatmap of mean gene expression per group.

    Aggregates mean expression of ``genes`` within each category of
    ``groupby`` and renders it as an interactive heatmap, optionally
    z-scored across rows or columns. Interactive mirror of
    :func:`complex_heatmap`.

    Parameters
    ----------
    adata : anndata.AnnData or str
        Annotated data matrix or path; aggregated via the internal
        ``_aggregate`` helper.
    genes : sequence of str
        Genes to include as heatmap columns.
    groupby : str
        Column in ``adata.obs`` whose categories form the heatmap rows.
    layer : str, optional
        Name of an ``adata.layers`` entry to aggregate. When ``None``
        (default) the main ``X`` matrix (or ``raw`` if ``use_raw``) is
        used.
    use_raw : bool, default False
        If ``True`` aggregate from ``adata.raw`` instead of ``adata.X``.
    z_score : str or None, default "row"
        Standardisation applied to the mean matrix: ``"row"`` z-scores
        each gene across groups, ``"col"`` z-scores each group across
        genes, and ``None`` leaves raw means. When z-scoring the
        colourbar is centred at zero.
    cmap : str, default "RdBu_r"
        Plotly continuous colourscale name.
    title : str, optional
        Figure title. When ``None`` no title is shown.
    save : str, optional
        Path to save the figure to via ``_maybe_save``. When ``None`` the
        figure is not written to disk.
    height : int, default 500
        Minimum figure height in pixels; grows with the number of groups.
    width : int, default 800
        Minimum figure width in pixels; grows with the number of genes.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive heatmap figure.
    """
    import plotly.graph_objects as go
    from .complex_heatmap import _aggregate

    mean_df, _ = _aggregate(adata, groupby, list(genes), layer=layer, use_raw=use_raw)
    if z_score == "row":
        mean_df = (
            mean_df.sub(mean_df.mean(axis=1), axis=0)
            .div(mean_df.std(axis=1).replace(0, np.nan), axis=0)
            .fillna(0)
        )
    elif z_score == "col":
        mean_df = (
            mean_df.sub(mean_df.mean(axis=0), axis=1)
            .div(mean_df.std(axis=0).replace(0, np.nan), axis=1)
            .fillna(0)
        )
    fig = go.Figure(
        go.Heatmap(
            z=mean_df.values,
            x=list(mean_df.columns),
            y=list(mean_df.index),
            colorscale=cmap,
            zmid=0 if z_score in ("row", "col") else None,
            colorbar=dict(title="z-score" if z_score else "mean"),
            hovertemplate=f"{groupby}: %{{y}}<br>gene: %{{x}}<br>value: %{{z:.2f}}<extra></extra>",
        )
    )
    fig.update_layout(
        title=title,
        height=max(height, 150 + 22 * len(mean_df.index)),
        width=max(width, 200 + 26 * len(mean_df.columns)),
        xaxis_title="gene",
        yaxis_title=groupby,
        template="plotly_white",
        margin=dict(l=100, r=80, t=50, b=110),
    )
    fig.update_xaxes(tickangle=-45, automargin=True)
    fig.update_yaxes(automargin=True)
    _maybe_save(fig, save)
    return fig


# ---------------------------------------------------------------------------
# Complex dotplot (Plotly dot heatmap, size = % expressing, colour = mean)
# ---------------------------------------------------------------------------


def interactive_complex_dotplot(
    adata,
    genes: Sequence[str],
    groupby: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    expression_cutoff: float = 0,
    cmap: str = "Reds",
    title: Optional[str] = None,
    save: Optional[str] = None,
    height: int = 500,
    width: int = 800,
    size_max: int = 22,
):
    """Interactive Plotly dot heatmap of expression (size = % expressing, colour = mean).

    Aggregates, per category of ``groupby``, the mean expression and the
    fraction of expressing cells for each gene, then draws a grid of
    dots where dot size encodes the expressing fraction and dot colour
    encodes the mean expression. Interactive mirror of
    :func:`complex_dotplot`.

    Parameters
    ----------
    adata : anndata.AnnData or str
        Annotated data matrix or path; aggregated via the internal
        ``_aggregate`` helper.
    genes : sequence of str
        Genes to include as columns of the dot grid.
    groupby : str
        Column in ``adata.obs`` whose categories form the rows.
    layer : str, optional
        Name of an ``adata.layers`` entry to aggregate. When ``None``
        (default) the main ``X`` matrix (or ``raw`` if ``use_raw``) is
        used.
    use_raw : bool, default False
        If ``True`` aggregate from ``adata.raw`` instead of ``adata.X``.
    expression_cutoff : float, default 0
        Threshold above which a cell is counted as expressing a gene when
        computing the expressing-fraction (dot size).
    cmap : str, default "Reds"
        Plotly continuous colourscale name used for the mean-expression
        colour.
    title : str, optional
        Figure title. When ``None`` no title is shown.
    save : str, optional
        Path to save the figure to via ``_maybe_save``. When ``None`` the
        figure is not written to disk.
    height : int, default 500
        Minimum figure height in pixels; grows with the number of groups.
    width : int, default 800
        Minimum figure width in pixels; grows with the number of genes.
    size_max : int, default 22
        Maximum dot diameter in pixels, used to scale the expressing
        fraction (0-1) into marker sizes.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive dot-heatmap figure.
    """
    import plotly.graph_objects as go
    from .complex_heatmap import _aggregate

    mean_df, frac_df = _aggregate(
        adata,
        groupby,
        list(genes),
        layer=layer,
        use_raw=use_raw,
        expression_cutoff=expression_cutoff,
    )
    cats = list(mean_df.index)
    genes = list(mean_df.columns)
    yy, xx = np.meshgrid(np.arange(len(cats)), np.arange(len(genes)), indexing="ij")
    sizes = frac_df.values
    fig = go.Figure(
        go.Scatter(
            x=xx.ravel(),
            y=yy.ravel(),
            mode="markers",
            marker=dict(
                size=sizes.ravel() * size_max,
                sizemode="diameter",
                sizemin=1,
                color=mean_df.values.ravel(),
                colorscale=cmap,
                line=dict(width=0.5, color="black"),
                colorbar=dict(title="mean"),
            ),
            customdata=np.stack([mean_df.values.ravel(), sizes.ravel()], axis=-1),
            hovertemplate=(
                f"{groupby}: %{{text}}<br>gene: %{{meta}}<br>"
                "mean: %{customdata[0]:.2f}<br>fraction: %{customdata[1]:.2f}<extra></extra>"
            ),
            text=[cats[i] for i in yy.ravel()],
            meta=[genes[j] for j in xx.ravel()],
        )
    )
    fig.update_layout(
        title=title,
        height=max(height, 150 + 24 * len(cats)),
        width=max(width, 200 + 26 * len(genes)),
        template="plotly_white",
        xaxis=dict(
            tickmode="array",
            tickvals=np.arange(len(genes)),
            ticktext=genes,
            tickangle=-45,
            automargin=True,
        ),
        yaxis=dict(
            tickmode="array",
            tickvals=np.arange(len(cats)),
            ticktext=cats,
            autorange="reversed",
            automargin=True,
        ),
        margin=dict(l=120, r=80, t=50, b=110),
    )
    _maybe_save(fig, save)
    return fig


__all__ += [
    "interactive_dot_plot",
    "interactive_chord",
    "interactive_upset",
    "interactive_complex_heatmap",
    "interactive_complex_dotplot",
]


# ---------------------------------------------------------------------------
# Ridgeline (joyplot)
# ---------------------------------------------------------------------------


def _gene_values_by_group(adata, gene, groupby, layer=None, use_raw=False, query=None):
    """Return an ordered ``{category: np.ndarray}`` of a gene's values.

    Reads the single gene from a (possibly backed) AnnData, handling
    ``layer`` / ``use_raw`` selection and an optional ``obs.query`` subset,
    and groups the per-cell values by ``groupby``.
    """
    import anndata

    if isinstance(adata, str):
        adata = anndata.read_h5ad(adata, backed="r")
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("interactive_ridgeline requires an AnnData or .h5ad path.")

    rows = adata.obs.query(query).index if query is not None else None
    use_raw = use_raw and adata.raw is not None
    src_names = list(adata.raw.var_names if use_raw else adata.var_names)
    if gene not in src_names:
        raise KeyError(f"gene {gene!r} not in var_names.")

    # Read only the requested rows/gene so backed AnnData does not
    # materialise the full matrix. ``layer`` / ``use_raw`` need the gene
    # column resolved against the right var space.
    row_sel = rows if rows is not None else slice(None)
    if use_raw:
        # Raw lives on the parent object; subset rows first, then read raw.
        sub = adata[row_sel].to_memory() if adata.isbacked else adata[row_sel]
        col = list(sub.raw.var_names).index(gene)
        X = sub.raw.X[:, col]
        groups = sub.obs[groupby].astype(str).to_numpy()
    else:
        if adata.isbacked:
            sub = adata[row_sel].to_memory()[:, gene]
        else:
            sub = adata[row_sel, gene]
        if layer is not None and layer in sub.layers:
            X = sub.layers[layer][:, 0]
        else:
            X = sub.X[:, 0]
        groups = sub.obs[groupby].astype(str).to_numpy()
    if hasattr(X, "toarray"):
        X = X.toarray()
    values = np.asarray(X).ravel().astype(float)
    if adata.isbacked:
        try:
            adata.file.close()
        except Exception:
            pass
    out = {}
    for g in pd.unique(groups):
        out[g] = values[groups == g]
    return out


def interactive_ridgeline(
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
    fill_alpha: float = 0.85,
    linecolor: str = "white",
    linewidth: float = 1.0,
    title: Optional[str] = None,
    xaxis_title: Optional[str] = None,
    yaxis_title: Optional[str] = None,
    save: Optional[str] = None,
    height: Optional[int] = None,
    width: int = 420,
    query: Optional[str] = None,
):
    """Interactive ridgeline (joyplot) of a gene's expression per group.

    Plotly mirror of :func:`adataviz.pl.ridgeline`. Draws one smoothed KDE
    curve per category of ``groupby`` as a filled area, stacked vertically
    with a controllable overlap, so the distribution of ``gene`` can be
    compared across groups and inspected interactively (hover shows the
    group, expression value and density). Returns a
    :class:`plotly.graph_objects.Figure`.

    Parameters
    ----------
    adata : anndata.AnnData or str
        Annotated data matrix, or path to an ``.h5ad`` file (opened
        backed). Expression for ``gene`` is read from ``X`` /
        ``raw.X`` / ``layers[layer]``.
    gene : str
        Single gene to plot; must exist in the resolved var space.
    groupby : str
        Column in ``adata.obs`` used as the categorical grouping; each
        category becomes one ridge (row).
    layer : str, optional
        ``adata.layers`` entry to read expression from. Takes precedence
        over ``X`` when present.
    use_raw : bool, default False
        Read expression / gene names from ``adata.raw`` when available.
        Ignored when ``layer`` is given.
    order : sequence, optional
        Explicit ordering of ``groupby`` categories from bottom to top.
        Categories not listed are dropped. None uses the natural order.
    palette : dict, str or None, optional
        Colour mapping for the categories. Accepts a ``{category: colour}``
        mapping, an ``.xlsx`` palette path, or ``None`` to fall back to
        colours from ``adata.uns`` or a default.
    overlap : float, default 0.7
        Fraction of vertical overlap between adjacent ridges. ``0`` draws
        non-overlapping curves; larger values push each curve further into
        its neighbours (peak height equals ``1 + overlap`` rows).
    bw_method : str, float or None, optional
        Bandwidth selector forwarded to
        :class:`scipy.stats.gaussian_kde`.
    n_points : int, default 256
        Number of points on the shared x-grid where each density is
        evaluated.
    clip : tuple of float, optional
        ``(low, high)`` limits for the x-grid. None spans the observed
        min/max (with a small margin).
    fill_alpha : float, default 0.85
        Opacity of the filled area under each curve.
    linecolor : str, default "white"
        Colour of the density outline.
    linewidth : float, default 1.0
        Width of the density outline.
    title : str, optional
        Figure title. None uses ``"{gene} Expression"``.
    xaxis_title : str, optional
        X axis label. Defaults to ``"Normalized expression"``.
    yaxis_title : str, optional
        Y axis label. Defaults to ``groupby``.
    save : str, optional
        Output path. ``.html`` writes standalone HTML; any other extension
        writes a static image. None disables saving.
    height : int, optional
        Figure height in pixels. None derives it from the group count.
    width : int, default 420
        Figure width in pixels.
    query : str, optional
        Pandas-style query string applied to ``adata.obs`` to subset cells
        before plotting (efficient for backed AnnData).

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive ridgeline figure.
    """
    import plotly.graph_objects as go
    from scipy.stats import gaussian_kde

    by_group = _gene_values_by_group(
        adata, gene, groupby, layer=layer, use_raw=use_raw, query=query
    )
    present = pd.Series(list(by_group.keys()))
    cats = categorical_order(present, order)
    cats = [c for c in cats if c in by_group]

    _obs, ad = resolve_adata_obs(adata) if not isinstance(adata, str) else (None, None)
    colors = resolve_palette(
        palette, cats, sheet_name=groupby, adata=ad, groupby=groupby
    )

    all_vals = np.concatenate([by_group[c] for c in cats]) if cats else np.array([0.0])
    all_vals = all_vals[np.isfinite(all_vals)]
    if clip is None:
        lo, hi = float(np.min(all_vals)), float(np.max(all_vals))
        margin = 0.05 * (hi - lo or 1.0)
        clip = (lo - margin, hi + margin)
    grid = np.linspace(clip[0], clip[1], n_points)

    def _density(v):
        v = np.asarray(v, dtype=float)
        v = v[np.isfinite(v)]
        if v.size < 2 or np.allclose(v, v[0]):
            hist, edges = np.histogram(
                v, bins=max(10, n_points // 8), density=True
            )
            centers = 0.5 * (edges[:-1] + edges[1:])
            return np.interp(grid, centers, hist, left=0.0, right=0.0)
        return gaussian_kde(v, bw_method=bw_method)(grid)

    densities = {c: _density(by_group[c]) for c in cats}
    dmax = max((float(np.max(d)) for d in densities.values() if d.size), default=1.0)
    step = 1.0
    amp = (1.0 + max(0.0, overlap)) * step
    scale = amp / dmax if dmax > 0 else 1.0

    def _hex_to_rgba(color: str, alpha: float) -> str:
        if not isinstance(color, str) or not color.startswith("#"):
            return color
        h = color.lstrip("#")
        r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
        return f"rgba({r},{g},{b},{alpha})"

    n = len(cats)
    fig = go.Figure()
    # Draw top-to-bottom so lower ridges render in front of the ones above.
    for i, c in enumerate(reversed(cats)):
        idx = n - 1 - i  # actual row index (0 = bottom)
        base = idx * step
        d = densities[c]
        y = base + d * scale
        fig.add_trace(
            go.Scatter(
                x=np.concatenate([grid, grid[::-1]]),
                y=np.concatenate([y, np.full_like(grid, base)]),
                fill="toself",
                mode="lines",
                line=dict(color=linecolor, width=linewidth),
                fillcolor=_hex_to_rgba(colors[c], fill_alpha),
                name=str(c),
                hovertemplate=(
                    f"{groupby}: {c}<br>value: %{{x:.3g}}"
                    "<br>density: %{customdata:.3g}<extra></extra>"
                ),
                customdata=np.concatenate([d, d[::-1]]),
                showlegend=False,
            )
        )

    fig.update_layout(
        title=title if title is not None else f"{gene} Expression",
        template="plotly_white",
        width=width,
        height=height if height is not None else int(max(320, 42 * n + 120)),
        margin=dict(l=70, r=30, t=50, b=50),
    )
    fig.update_xaxes(
        title=xaxis_title if xaxis_title is not None else "Normalized expression",
        range=[clip[0], clip[1]],
        zeroline=False,
    )
    fig.update_yaxes(
        title=yaxis_title if yaxis_title is not None else groupby,
        tickmode="array",
        tickvals=[i * step for i in range(n)],
        ticktext=cats,
        range=[-0.5 * step, (n - 1) * step + amp + 0.3 * step],
        showgrid=False,
        zeroline=False,
    )
    _maybe_save(fig, save)
    return fig


__all__ += ["interactive_ridgeline"]
