#!/usr/bin/env python
# -*- coding: utf-8 -*-
import importlib as _importlib
import os
import fire
from ._version import version as __version__  # noqa: F401


def __getattr__(name):
    """
    Lazy module loader for adataviz submodules.

    Loads ``tl`` (tools) and ``pl`` (plotting) submodules on first access
    to avoid importing heavy dependencies (plotly, scanpy, seaborn)
    at package import time.

    Parameters
    ----------
    name : str
            Attribute name. Supports "tl" for tools module and "pl" for
            plotting module.

    Returns
    -------
    module
            The requested submodule.

    Raises
    ------
    AttributeError
            If ``name`` is not "tl" or "pl".
    """
    if name == "tl":
        return _importlib.import_module(".tools", __name__)
    if name == "pl":
        return _importlib.import_module(".plotting", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


_ROOT = os.path.abspath(os.path.dirname(__file__))


def adataviz(command=None):
    """
    Usage: adataviz command [options]\n

    Available commands are:
            plot:	plot adata
            tool:	Other tools

    Parameters
    ----------
    command :
            choice from ['plot', 'tool']

    Returns
    -------

    """
    doc_string = """
Available subcommands:
	plot:	plot adata
	tool:	Other tools
		"""
    if command is None:
        return doc_string
    command = command.lower()
    _pl = _importlib.import_module(".plotting", __name__)
    _tl = _importlib.import_module(".tools", __name__)
    if command == "plot":
        return {
            "plot_categorical": _pl.plot_categorical,
            "plot_genes": _pl.plot_genes,
        }
    elif command == "tool":
        return {
            "scrna2pseudobulk": _tl.scrna2pseudobulk,
            "stat_pseudobulk": _tl.stat_pseudobulk,
            "normalize_adata": _tl.normalize_adata,
            "export_pseudobulk_adata": _tl.export_pseudobulk_adata,
            "parse_gtf": _tl.parse_gtf,
            "downsample_adata": _tl.downsample_adata,
            "get_obs": _tl.get_obs,
            "composition": _tl.composition,
            "taxonomy": _tl.taxonomy,
            "get_markers": _tl.get_markers,
        }
    else:
        print(doc_string)
        exit()


def main():
    """CLI entry point. Dispatches commands via Python Fire."""
    from .utils import serialize

    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire(adataviz, serialize=serialize)


if __name__ == "__main__":
    main()
