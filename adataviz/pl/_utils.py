"""Internal helpers shared across :mod:`adataviz.pl` modules.

These utilities normalise inputs (AnnData / DataFrame / file paths),
resolve colour palettes, and handle figure save/show in a way that is
forward-compatible with new pandas / matplotlib / anndata releases (no
deprecated ``is_categorical_dtype``, no removed positional kwargs).
"""

from __future__ import annotations

import contextlib
import os
import warnings
from typing import Any, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

__all__ = [
    "resolve_obs",
    "resolve_adata_obs",
    "resolve_palette",
    "to_hex",
    "save_or_show",
    "categorical_order",
    "categorical_series",
    "gene_positions",
    "merge_legend_kws",
    "maybe_adjust_texts",
    "palette_from_uns",
    "use_scientific_style",
    "show_fig",
    "import_pch",
    "quiet_pch",
    "apply_layout",
    "label_extent",
    "outside_legend",
]


def gene_positions(var_names, genes: Sequence[str]) -> list:
    """Column positions of ``genes`` in ``var_names`` (``KeyError`` if absent).

    Uses the fast ``Index.get_indexer`` path and falls back to a dict lookup
    (last occurrence wins, matching the previous behaviour) when
    ``var_names`` is not unique, which ``get_indexer`` refuses.
    """
    var = pd.Index(var_names)
    if var.is_unique:
        idx = var.get_indexer(list(genes))
    else:
        pos = {g: i for i, g in enumerate(var)}
        idx = np.array([pos.get(g, -1) for g in genes], dtype=int)
    if (idx < 0).any():
        missing = [g for g, i in zip(genes, idx) if i < 0]
        raise KeyError(f"genes not in adata.var_names: {missing[:5]}")
    return idx.tolist()


def import_pch():
    """Import :mod:`PyComplexHeatmap` with pandas>=3 compatibility shims.

    pandas 3 removed ``DataFrame.applymap`` (renamed ``DataFrame.map`` in
    2.1), but PyComplexHeatmap <= 1.8.x still calls it inside
    ``DotClustermapPlotter``. Alias it back so the dot plots keep working
    on new pandas without pinning the upstream package.
    """
    if not hasattr(pd.DataFrame, "applymap") and hasattr(pd.DataFrame, "map"):
        pd.DataFrame.applymap = pd.DataFrame.map  # type: ignore[attr-defined]
    import PyComplexHeatmap as pch

    return pch


@contextlib.contextmanager
def quiet_pch(fig):
    """Silence upstream PyComplexHeatmap warnings while it draws ``fig``.

    - Disables any auto layout engine on the target figure so that
      PyComplexHeatmap's manual GridSpec layout does not trigger the
      "Axes are not compatible with tight_layout" warning on canvas draw.
    - Filters the pandas ``FutureWarning`` about callable ``mean`` that
      comes from ``PyComplexHeatmap.dotHeatmap``.
    """
    # Force this figure off any auto-layout engine (tight/constrained) a user
    # rcParams ``figure.autolayout=True`` would otherwise install: PCH lays
    # its GridSpec out by hand.
    try:
        fig.set_layout_engine("none")
    except Exception:
        pass
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=FutureWarning,
            module=r"PyComplexHeatmap\..*",
        )
        warnings.filterwarnings(
            "ignore",
            message=r".*not compatible with tight_layout.*",
            category=UserWarning,
        )
        yield


def apply_layout(fig, own: bool = True, **engine_kws):
    """Make ``fig`` lay itself out at draw time with the constrained engine.

    Constrained layout (unlike ``tight_layout``) accounts for legends placed
    outside their axes via ``bbox_to_anchor`` as well as rotated tick labels
    and figure titles, so panels never overlap their decorations. Only
    applied when ``own`` is true (the figure was created by us); a figure the
    caller handed in through ``ax=`` is left untouched. Extra keyword
    arguments (``h_pad``, ``w_pad``, ``hspace``, ``wspace``) are forwarded to
    the engine. Falls back to a one-off ``tight_layout`` on very old
    matplotlib.
    """
    if not own:
        return
    try:
        fig.set_layout_engine("constrained")
        if engine_kws:
            fig.get_layout_engine().set(**engine_kws)
    except Exception:  # pragma: no cover - matplotlib < 3.6
        try:
            fig.tight_layout()
        except Exception:
            pass


def label_extent(labels: Sequence, fontsize: float = 8, rotation: float = 0) -> float:
    """Rough inch-extent of the longest label (for auto figure sizing).

    Uses an average glyph width of ~0.6 em; ``rotation`` (degrees) projects
    the text onto the horizontal axis so a -45 degree tick label reserves
    about 70 % of its unrotated width.
    """
    import math

    maxlen = max((len(str(x)) for x in labels), default=0)
    width_in = maxlen * 0.6 * fontsize / 72.0
    return width_in * abs(math.cos(math.radians(rotation))) + 0.1


def outside_legend(fig, ax, handles=None, labels=None, legend_kws=None, **defaults):
    """Draw a legend that constrained layout keeps clear of *all* axes content.

    When the caller did not pin a location (``loc`` / ``bbox_to_anchor`` in
    ``legend_kws``) and matplotlib supports it, the legend is attached to the
    figure with ``loc="outside center right"`` so it never collides with
    labels that stick out of the axes (polar tick labels, rotated x labels,
    ...). Otherwise it is a regular ``ax.legend`` with the usual
    right-of-axes defaults.
    """
    kws = merge_legend_kws(legend_kws, **defaults)
    user_pinned = bool(legend_kws) and (
        "loc" in legend_kws or "bbox_to_anchor" in legend_kws
    )
    if handles is not None:
        kws["handles"] = handles
    if labels is not None:
        kws["labels"] = labels
    if not user_pinned and fig.get_layout_engine() is not None:
        kws.pop("bbox_to_anchor", None)
        kws["loc"] = "outside center right"
        try:
            return fig.legend(**kws)
        except (ValueError, TypeError):  # matplotlib < 3.7: no "outside" locs
            kws["loc"] = "center left"
            kws["bbox_to_anchor"] = (1.02, 0.5)
    return ax.legend(**kws)


def resolve_obs(obs: Any) -> pd.DataFrame:
    """Return a :class:`pandas.DataFrame` from an obs-like input.

    Accepts a DataFrame, an :class:`anndata.AnnData`, or a path to a
    ``.h5ad`` / ``.csv`` / ``.tsv`` / ``.tsv.gz`` file. Always returns a
    copy.
    """
    if isinstance(obs, pd.DataFrame):
        return obs.copy()
    try:
        import anndata as _ad

        if isinstance(obs, _ad.AnnData):
            return obs.obs.copy()
    except ImportError:  # pragma: no cover
        pass
    if not isinstance(obs, (str, os.PathLike)):
        raise TypeError(
            f"obs must be a DataFrame, AnnData, or path; got {type(obs).__name__}"
        )
    path = os.path.abspath(os.path.expanduser(str(obs)))
    if path.endswith(".h5ad"):
        import anndata as _ad

        ad = _ad.read_h5ad(path, backed="r")
        try:
            df = ad.obs.copy()
        finally:
            try:
                ad.file.close()
            except Exception:
                pass
        return df
    if path.endswith(".csv"):
        return pd.read_csv(path, index_col=0)
    return pd.read_csv(path, sep="\t", index_col=0)


def resolve_adata_obs(adata: Any):
    """Return ``(obs_df, adata_or_none)`` for input that may be AnnData.

    Plot functions use this so they can both (a) operate on ``adata.obs``
    and (b) reach into ``adata.uns`` for colour palettes when available.
    Falls back to ``(_resolve_obs(adata), None)`` otherwise.
    """
    try:
        import anndata as _ad
    except ImportError:  # pragma: no cover
        return resolve_obs(adata), None
    if isinstance(adata, _ad.AnnData):
        return adata.obs.copy(), adata
    return resolve_obs(adata), None


def palette_from_uns(adata, groupby: str, categories: Sequence[str]) -> dict:
    """Build a ``{cat: hex}`` dict from ``adata.uns[f"{groupby}_colors"]``.

    Returns an empty dict if the entry is missing or malformed.
    """
    if adata is None:
        return {}
    key = f"{groupby}_colors"
    if key not in adata.uns or groupby not in adata.obs:
        return {}
    cats = list(adata.obs[groupby].cat.categories) if hasattr(
        adata.obs[groupby], "cat"
    ) else list(pd.unique(adata.obs[groupby]))
    colors = list(adata.uns[key])
    out = {str(c): col for c, col in zip(cats, colors)}
    # Filter to only requested categories.
    return {c: out[str(c)] for c in categories if str(c) in out}


def resolve_palette(
    palette: Union[None, Mapping[str, str], str],
    categories: Sequence[str],
    *,
    sheet_name: Optional[str] = None,
    adata=None,
    groupby: Optional[str] = None,
) -> dict:
    """Return a ``{category: hex}`` dict from a flexible palette spec.

    Resolution order:

    1. Explicit dict / Excel-path entries.
    2. ``adata.uns[f"{groupby}_colors"]`` if *adata* and *groupby* given.
    3. matplotlib ``tab10`` / ``tab20`` defaults for any remaining
       categories.
    """
    cats = list(dict.fromkeys(str(c) for c in categories))
    base: dict = {}
    if isinstance(palette, Mapping):
        base = {str(k): v for k, v in palette.items()}
    elif isinstance(palette, str):
        path = os.path.abspath(os.path.expanduser(palette))
        if os.path.exists(path):
            sheet = sheet_name if sheet_name is not None else 0
            try:
                df = pd.read_excel(path, sheet_name=sheet, index_col=0)
                base = df["Hex"].to_dict()
            except Exception:
                base = {}
    if adata is not None and groupby is not None:
        for c, col in palette_from_uns(adata, groupby, cats).items():
            base.setdefault(c, col)
    cmap_name = "tab20" if len(cats) > 10 else "tab10"
    cmap = plt.get_cmap(cmap_name)
    out: dict = {}
    for i, c in enumerate(cats):
        out[c] = base.get(c, base.get(str(c), to_hex(cmap(i % cmap.N))))
    return out


def to_hex(rgba) -> str:
    from matplotlib.colors import to_hex as _to_hex

    return _to_hex(rgba)


def save_or_show(fig, save: Optional[str], show: bool = False):
    """Save a Matplotlib :class:`Figure` and optionally display it."""
    if save:
        out = os.path.abspath(os.path.expanduser(save))
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        fig.savefig(out, bbox_inches="tight")
    if show:
        plt.show()


def categorical_order(values: pd.Series, order: Optional[Sequence] = None) -> list:
    """Return ordered list of categories in *values*, dropping NaN.

    When ``order`` is provided, returns that order (filtered to values
    that occur in *values*). Otherwise uses categorical order or the
    sorted unique values.
    """
    ser = values.dropna()
    present = set(map(str, ser.unique()))
    if order is not None:
        return [str(o) for o in order if str(o) in present]
    if isinstance(ser.dtype, pd.CategoricalDtype):
        return [str(c) for c in ser.cat.categories if str(c) in present]
    return sorted(present)


def categorical_series(values: pd.Series, order: Optional[Sequence] = None) -> pd.Series:
    """Return *values* as an ordered ``Categorical`` of ``str`` labels.

    Keeps the categorical order of the input (or ``order`` when given, or the
    sorted unique values for plain object/numeric columns) while making the
    labels plain strings, so downstream ``groupby`` results come back in a
    meaningful order instead of alphabetically. Rows whose category is not in
    the final order become ``NaN`` and are dropped by ``observed=True``
    group-bys.
    """
    cats = categorical_order(values, order)
    return pd.Series(
        pd.Categorical(values.astype(str), categories=cats, ordered=True),
        index=values.index,
        name=values.name,
    )


_DEFAULT_LEGEND_KWS = dict(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    fontsize=8,
    handlelength=1.2,
    handletextpad=0.5,
    borderpad=0.4,
    labelspacing=0.4,
)


def merge_legend_kws(user_kws: Optional[Mapping[str, Any]], **defaults) -> dict:
    """Merge user-supplied legend kwargs with module defaults.

    Defaults follow ``plot_genes`` (compact, right of the axes). The user
    dict wins for any key it provides.
    """
    out = dict(_DEFAULT_LEGEND_KWS)
    out.update(defaults)
    if user_kws:
        out.update(user_kws)
    return out


def maybe_adjust_texts(texts, ax=None, **kwargs):
    """Optionally call :func:`adjustText.adjust_text` if installed.

    Returns ``True`` if adjustment was performed, ``False`` if the
    optional dependency is missing.
    """
    if not texts:
        return False
    try:
        from adjustText import adjust_text  # type: ignore
    except ImportError:
        return False
    kw = dict(
        only_move={"text": "y", "static": "y", "explode": "y", "pull": "y"},
        arrowprops=dict(arrowstyle="-", color="grey", lw=0.5, alpha=0.6),
        ax=ax,
    )
    kw.update(kwargs)
    if kw.get("ax") is None:
        kw.pop("ax", None)
    adjust_text(texts, **kw)
    return True


# ---------------------------------------------------------------------------
# Style + plotly helpers
# ---------------------------------------------------------------------------


def use_scientific_style():
    """Configure matplotlib rcParams for publication-quality figures."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": "Arial",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "figure.titlesize": 11,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.title_fontsize": 9,
            "lines.linewidth": 1.2,
            "axes.linewidth": 1.0,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "figure.figsize": (6.5, 4.5),
            # NOTE: ``figure.autolayout`` is deliberately *not* enabled: a
            # global tight_layout at draw time fights PyComplexHeatmap's
            # manual GridSpec and squeezes axes under legends anchored outside
            # them. Every adataviz plot manages its own layout instead.
            "figure.autolayout": False,
            "savefig.transparent": True,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def show_fig(fig, filename: str = "plot"):
    """Display a Plotly figure with an editable interactive toolbar.

    Renders ``fig`` with a customised Plotly mode-bar configuration that
    enables interactive editing (dragging titles/legends/annotations,
    scroll-to-zoom, etc.) and sets SVG as the default export format.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to display.
    filename : str, default ``"plot"``
        Base filename used when the figure is exported via the mode-bar
        "download" (camera) button; the image is saved as ``<filename>.svg``.

    Returns
    -------
    None
        The figure is shown in-place via ``fig.show(config=...)``.
    """
    config = {
        "displayModeBar": "hover",
        "showLink": False,
        "linkText": "Edit on plotly",
        "scrollZoom": True,
        "displaylogo": False,
        "toImageButtonOptions": {"format": "svg", "filename": filename},
        "modeBarButtonsToRemove": ["sendDataToCloud"],
        "editable": True,
        "autosizable": True,
        "edits": {
            "titleText": True,
            "legendPosition": True,
            "colorbarTitleText": True,
            "shapePosition": True,
            "annotationPosition": True,
            "annotationText": True,
            "axisTitleText": True,
            "legendText": True,
            "colorbarPosition": True,
        },
    }
    fig.show(config=config)
