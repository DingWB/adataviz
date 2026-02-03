#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
# from .clustermap import (heatmap, ClusterMapPlotter, composite,
#                          DendrogramPlotter)
# from .oncoPrint import oncoprint, oncoPrintPlotter
# __all__=['*']
import sys
from ._version import version as __version__
from . import tools as tl
from . import plotting as pl
import fire

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
	doc_string="""
Available subcommands:
	plot:	plot adata
	tool:	Other tools
		"""
	if command is None:
		return doc_string
	command=command.lower()
	if command == "plot":
		return {
			'plot_cluster': pl.plot_cluster,
			'plot_gene': pl.plot_gene,
			'plot_genes': pl.plot_genes,
		}
	elif command=="tool":
		return {
			"scrna2pseudobulk": tl.scrna2pseudobulk,
			'stat_pseudobulk':tl.stat_pseudobulk,
			'export_pseudobulk_adata': tl.export_pseudobulk_adata,
			'parse_gtf':tl.parse_gtf,
			'downsample_adata':tl.downsample_adata,
			'get_obs':tl.get_obs,
		}
	else:
		print(doc_string)
		exit()
		
def main():
	from .utils import serialize
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(adataviz, serialize=serialize)

if __name__=="__main__":
	main()