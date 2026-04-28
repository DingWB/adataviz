"""Internal helpers shared across :mod:`adataviz.pl` modules.

These utilities normalise inputs (AnnData / DataFrame / file paths),
resolve colour palettes, and handle figure save/show in a way that is
forward-compatible with new pandas / matplotlib / anndata releases (no
deprecated ``is_categorical_dtype``, no removed positional kwargs).
"""

from __future__ import annotations

import os
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
    "merge_legend_kws",
    "maybe_adjust_texts",
    "palette_from_uns",
    "use_scientific_style",
    "show_fig",
]


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
    if key not in adata.uns:
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
            "figure.autolayout": True,
            "savefig.transparent": True,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def show_fig(fig, filename: str = "plot"):
    """Display a Plotly figure with an editable interactive toolbar."""
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
