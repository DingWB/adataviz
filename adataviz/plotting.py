import os, sys
import pandas as pd
import anndata
import fire
import matplotlib.pylab as plt
import numpy as np
from matplotlib.colors import Normalize
import seaborn as sns
from .utils import (
	_make_tiny_axis_label, despine,
	zoom_ax, _extract_coords,
	_density_based_sample,_auto_size,
	_take_data_series, level_one_palette, 
	tight_hue_range,_text_anno_scatter,
	density_contour,plot_color_dict_legend,
	plot_marker_legend,plot_text_legend,plot_cmap_legend
)
from PyComplexHeatmap import HeatmapAnnotation,anno_label,anno_simple,DotClustermapPlotter
from .utils import normalize_mc_by_cell

def stacked_barplot(Input="cell_metadata_with_annotation.csv",groupby='Age',
					column='CellClass',x_order=None,y_order=None,linewidth=0.1,
					palette="~/Projects/mouse_pfc/metadata/mpfc_color_palette.xlsx",
					width=None,height=None,xticklabels_kws=None,outdir="figures",
					outname=None,lgd_kws=None,gap=0.05,sort_by=None):
	"""
	Plot stacked barplto to show the cell type composition in each `groupby` (
		such as Age, brain regions and so on.)
		For example: stacked_barplot(column='MajorType',width=3.5,height=6)
							stacked_barplot(column='CellClass',width=3.5,height=3)

	Parameters
	----------
	Input :
	groupby :
	column :
	x_order :
	y_order :
	linewidth :
	palette :
	width :
	height :
	xticklabels_kws :
	outdir :
	lgd_kws: dict

	Returns
	-------

	"""
	outdir=os.path.abspath(os.path.expanduser(outdir))
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	if isinstance(Input,pd.DataFrame):
		data=Input.copy()
	elif isinstance(Input, str) and Input.endswith('.h5ad'):
		input_path = os.path.abspath(os.path.expanduser(Input))
		adata = anndata.read_h5ad(input_path,backed='r')
		data=adata.obs
		del adata
	elif Input.endswith('.csv'):
		input_path = os.path.abspath(os.path.expanduser(Input))
		data=pd.read_csv(input_path,index_col=0)
	else:
		input_path = os.path.abspath(os.path.expanduser(Input))
		data = pd.read_csv(input_path, index_col=0,sep='\t')
	if outname is None:
		outname=os.path.join(outdir,f"{groupby}.{column}.barplot.pdf")
	xticklabels_kws={} if xticklabels_kws is None else xticklabels_kws
	xticklabels_kws.setdefault('rotation',-45)
	xticklabels_kws.setdefault("rotation_mode", "anchor")
	xticklabels_kws.setdefault('horizontalalignment', 'left') #see ?matplotlib.axes.Axes.tick_params
	xticklabels_kws.setdefault('verticalalignment', 'center')
	if not palette is None:
		if isinstance(palette,dict):
			color_palette=palette.copy()
		elif isinstance(palette,str) and os.path.exists(os.path.expanduser(palette)):
			palette=os.path.abspath(os.path.expanduser(palette))
			D=pd.read_excel(palette,
							sheet_name=None, index_col=0)
			color_palette=D[column].Hex.to_dict()
			keys=list(color_palette.keys())
			for k in data[column].unique():
				if k not in keys:
					color_palette[k]='gray'
		else:
			color_palette = palette
	else:
		color_palette=None
	df=data.groupby(groupby)[column].value_counts(normalize=True).unstack(level=column)
	if not sort_by is None:
		df.sort_values(sort_by,ascending=True,inplace=True)
	else:
		if x_order is None:
			x_order = sorted(df.index.tolist())
		if y_order is None:
			y_order = sorted(df.columns.tolist())
		df=df.loc[x_order,y_order]
	if width is None:
		width=max(df.shape[0]*0.45,10)
		if width < 2.5:
			width=2.5
	if height is None:
		height = max(df.shape[1]*0.5, 8)
		if height < 3.5:
			height = 3.5
	plt.figure()
	ax=df.plot.bar(stacked=True,align='edge', width=1-gap,edgecolor='black',
				   linewidth=linewidth,figsize=(width,height),
				   color=color_palette)
	ax.set_xlim(0,df.shape[0])
	ax.set_ylim(0,1)
	labels=[tick.get_text() for tick in ax.get_xticklabels()]
	ax.set_xticks(ticks=np.arange(0.5,df.shape[0],1),
				  labels=labels,**xticklabels_kws) # ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0.5,df.shape[1],1))) #ticker.MultipleLocator(0.5)
	ax.xaxis.label.set_visible(False)
	ax.tick_params(
		axis="y", #both
		which="both",left=False,right=False,labelleft=False,labelright=False,
		top=False,labeltop=False,#bottom=False,labelbottom=False
		)
	# ax.xaxis.set_tick_params(axis='x')
	lgd_kws = lgd_kws if not lgd_kws is None else {}  # bbox_to_anchor=(x,-0.05)
	lgd_kws.setdefault("frameon", True)
	lgd_kws.setdefault("ncol", 1)
	lgd_kws["loc"] = "upper left"
	lgd_kws.setdefault("borderpad", 0.1 * (1 / 25.4) * 72)  # 0.1mm
	lgd_kws.setdefault("markerscale", 1)
	lgd_kws.setdefault("handleheight", 1)  # font size, units is points
	lgd_kws.setdefault("handlelength", 1)  # font size, units is points
	lgd_kws.setdefault("borderaxespad", 0.1)  # The pad between the axes and legend border, in font-size units.
	lgd_kws.setdefault("handletextpad", 0.3)  # The pad between the legend handle and text, in font-size units.
	lgd_kws.setdefault("labelspacing", 0.1)  # gap height between two Patches,  0.05*mm2inch*72
	lgd_kws.setdefault("columnspacing", 1)
	lgd_kws.setdefault("bbox_to_anchor", (1, 1))
	lgd_kws.setdefault("title",column)
	ax.legend(**lgd_kws)
	plt.savefig(outname, transparent=True,bbox_inches='tight',dpi=300)
	plt.show()

def pieplot(Input="cell_metadata_with_annotation.csv",groupby='Age',outdir="figures",
			palette_path="~/Projects/mouse_pfc/metadata/mpfc_color_palette.xlsx",
			order=None,explode=0.05):
	outdir = os.path.abspath(os.path.expanduser(outdir))
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	# colors=None
	if isinstance(Input,pd.DataFrame):
		data=Input.copy()
	elif isinstance(Input, str) and Input.endswith('.h5ad'):
		input_path = os.path.abspath(os.path.expanduser(Input))
		print(f"Reading adata: {Input}")
		adata = anndata.read_h5ad(input_path, backed='r')
		# if f'{groupby}_colors' in adata.uns:
		#     colors={k:v for k,v in zip(adata.obs[groupby].cat.categories.tolist(),
		#                                adata.uns[f'{groupby}_colors'])}
		# else:
		#     colors=None
		data = adata.obs
		del adata
	elif Input.endswith('.csv'):
		input_path = os.path.abspath(os.path.expanduser(Input))
		data = pd.read_csv(input_path, index_col=0)
	else:
		input_path = os.path.abspath(os.path.expanduser(Input))
		data = pd.read_csv(input_path, index_col=0, sep='\t')

	if not palette_path is None:
		palette_path=os.path.abspath(os.path.expanduser(palette_path))
		D=pd.read_excel(palette_path,
						sheet_name=None, index_col=0)
		color_palette=D[groupby].Hex.to_dict()
	else:
		color_palette=None
	D=data[groupby].value_counts()
	if order is None:
		order=list(sorted(D.keys()))

	plt.figure()
	plt.pie([D[k] for k in order], labels=order,
			colors=[color_palette[k] for k in order],
			explode=[explode]*len(order), autopct='%.1f%%')
	# Add title to the chart
	plt.title('Distribution of #of cells across different stages')
	plt.savefig(os.path.join(outdir, groupby + '.piechart.pdf'), transparent=True,bbox_inches='tight',dpi=300)
	plt.show()

def plot_pseudotime(
	pseudotime="dpt_pseudotime.tsv",groupby='Age',y='dpt_pseudotime',
	hue=None,figsize=(5,3.5),outdir="figures",rotate=None,ylabel='Pseudotime',
	palette_path="~/Projects/mouse_pfc/metadata/mpfc_color_palette.xlsx",
):
	"""
	Plot pseudotime. plot_pseudotime(figsize=(6,3.5),groupby='MajorType',
									rotate=-45);
								plot_pseudotime(figsize=(3.5,3),groupby='CellClass')
								plot_pseudotime(figsize=(3.5,3),groupby='Age')

	Parameters
	----------
	pseudotime :
	groupby :
	y :
	hue :
	figsize :
	outdir :
	rotate :
	palette_path :

	Returns
	-------

	"""
	outdir = os.path.abspath(os.path.expanduser(outdir))
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	if not palette_path is None:
		palette_path=os.path.abspath(os.path.expanduser(palette_path))
		D=pd.read_excel(palette_path,
						sheet_name=None, index_col=0)
		color_palette=D[groupby].Hex.to_dict()
	else:
		color_palette=None

	data=pd.read_csv(os.path.expanduser(pseudotime),sep='\t',index_col=0)
	data.dpt_pseudotime.replace({np.inf: 1},inplace=True)
	order=data.groupby(groupby)[y].mean().sort_values(ascending=True).index.tolist()
	if not hue is None:
		hue_order=data.groupby(hue)[y].mean().sort_values(ascending=True).index.tolist()
	else:
		hue_order=None
	plt.figure(figsize=figsize)
	# ax = sns.swarmplot(data=data, palette=color_palette, \
	# 				   edgecolor='white', x=groupby, y=y, \
	# 				   order=order)
	ax=sns.violinplot(data=data, x=groupby, y=y,
				   scale='count', bw=.2, inner=None, saturation=0.6,
				   palette=color_palette, order=order,hue=hue,hue_order=hue_order)
	# plt.legend(frameon=True)
	ax.set_ylabel(ylabel)
	if not rotate is None:
		plt.setp(ax.xaxis.get_majorticklabels(), rotation=rotate,
				 rotation_mode="anchor",horizontalalignment='left')
	outname=groupby + '.pseudotime_violin.pdf' if hue is None else groupby + f'_{hue}.pseudotime_violin.pdf'
	plt.savefig(os.path.join(outdir, outname), transparent=True,bbox_inches='tight',dpi=300)
	plt.show()

def stacked_violinplot(adata, use_genes=None, groupby='Age',
					   cell_groups=None, parent=None, figsize=(6, 4), cmap='viridis'):
	import scanpy as sc
	ax = sc.pl.stacked_violin(
            adata[adata.obs[cell_groups[0]]==parent],
            var_names=use_genes, title=use_key, colorbar_title="Avg mc frac",
            groupby=groupby, dendrogram=True, swap_axes=False,
            cmap=cmap,
            figsize=figsize, scale='count', standard_scale='obs', inner='quart',
            # stripplot=False,jitter=False,
            show=False, layer=None)
	ax1 = ax['mainplot_ax']
	ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
	ax1.yaxis.set_tick_params(which='minor',left=True)
	ax1.grid(axis='y', linestyle='--', color='black',
				alpha=1, zorder=-5, which='minor')
	plt.savefig(f"{fig_basename}.{groupby}.stacked_violin.pdf",
				transparent=True,bbox_inches='tight',dpi=300)
	plt.show()

def categorical_scatter(
    data,ax=None,
    coord_base="umap",x=None,y=None, # coords
    hue=None,palette="auto",color=None, # color
    text_anno=None,text_kws=None,luminance=None,text_transform=None,
    dodge_text=False,dodge_kws=None, # text annotation
    show_legend=False,legend_kws=None, # legend
    s="auto",size=None,sizes=None, # sizes is a dict
    size_norm=None,size_portion=0.95, 
    axis_format="tiny",max_points=50000,
    labelsize=4,linewidth=0.5,zoomxy=1.05,
    outline=None,outline_pad=3,alpha=0.7,
    outline_kws=None,scatter_kws=None,
    rasterized="auto",coding=False,
	id_marker=True,legend_color_text=True,
	rectangle_marker=False,marker_fontsize=4,marker_pad=0.1,
):
	"""
	This function was copied from ALLCools and made some modifications.
	Plot categorical scatter plot with versatile options.

	Parameters
	----------
	rasterized
		Whether to rasterize the figure.
	return_fig
		Whether to return the figure.
	size_portion
		The portion of the figure to be used for the size norm.
	data
		Dataframe that contains coordinates and categorical variables
	ax
		this function do not generate ax, must provide an ax
	coord_base
		coords name, if provided, will automatically search for x and y
	x
		x coord name
	y
		y coord name
	hue : str
		categorical col name or series for color hue.
	palette : str or dict
		palette for color hue.
	color
		specify single color for all the dots
	text_anno
		categorical col name or series for text annotation.
	text_kws
		kwargs pass to plt.text, see: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
		including bbox, to see parameter for bbox, go to: https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.FancyBboxPatch.html#matplotlib.patches.FancyBboxPatch
		commonly used parameters are: 
		```
		text_kws=dict(fontsize=5,fontweight='black',
					color='black', #color could be a dict, keys are text to be annotated
					bbox=dict(boxstyle='round',edgecolor=(0.5, 0.5, 0.5, 0.2),fill=False,
								facecolor=(0.8, 0.8, 0.8, 0.2), #facecolor could also be a dict
								alpha=1,linewidth=0.5)
					)
		```
	text_transform
		transform for text annotation.
	dodge_text
		whether to dodge text annotation.
	dodge_kws
		kwargs for dodge text annotation.
	show_legend
		whether to show legend.
	legend_kws
		kwargs for legend.
	s
		single size value of all the dots.
	size
		mappable size of the dots.
	sizes
		mapping size to the sizes value.
	size_norm
		normalize size range for mapping.
	axis_format
		axis format.
	max_points
		maximum number of points to plot.
	labelsize
		label size pass to `ax.text`
	linewidth
		line width pass to `ax.scatter`
	zoomxy
		zoom factor for x and y-axis.
	outline
		categorical col name or series for outline.
	outline_pad
		outline padding.
	outline_kws
		kwargs for outline.
	scatter_kws
		kwargs for scatter.

	Returns
	-------
	if return_fig is True, return the figure and axes.
	else, return None.
	"""
	if ax is None:
		# fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
		ax = plt.gca()
	# add coords
	_data, x, y = _extract_coords(data, coord_base, x, y)
	# _data has 2 cols: "x" and "y", index are obs_names

	# down sample plot data if needed.
	if max_points is not None:
		if _data.shape[0] > max_points:
			_data = _density_based_sample(_data, seed=1, size=max_points, coords=["x", "y"])
	n_dots = _data.shape[0]

	# determine rasterized
	if rasterized == "auto":
		if n_dots > 200:
			rasterized = True
		else:
			rasterized = False

	# auto size if user didn't provide one
	if s == "auto":
		s = _auto_size(ax, n_dots)

	# default scatter options
	_scatter_kws = {"linewidth": 0, "s": s, "legend": None, "palette": palette, "rasterized": rasterized}
	if color is not None:
		if hue is not None:
			raise ValueError("Only one of color and hue can be provided")
		_scatter_kws["color"] = color
	if scatter_kws is not None:
		_scatter_kws.update(scatter_kws)

	# deal with color
	palette_dict = None
	if hue is not None:
		if isinstance(hue, str):
			_data["hue"] = _take_data_series(data, hue)
		else:
			_data["hue"] = hue.copy()
		_data["hue"] = _data["hue"].astype("category").cat.remove_unused_categories()

		# if the object has get_palette method, use it (AnnotZarr)
		palette = _scatter_kws["palette"]
		# deal with other color palette
		if palette_dict is None:
			if isinstance(palette, str) or isinstance(palette, list):
				palette_dict = level_one_palette(_data["hue"], order=None, palette=palette)
			elif isinstance(palette, dict):
				palette_dict = palette
			else:
				raise TypeError(f"Palette can only be str, list or dict, " f"got {type(palette)}")
		_scatter_kws["palette"] = palette_dict

	# deal with size
	if size is not None:
		if isinstance(size, str):
			_data["size"] = _take_data_series(data, size).astype(float)
		else:
			_data["size"] = size.astype(float)
		size = "size"

		if size_norm is None:
			# get the smallest range that include "size_portion" of data
			size_norm = tight_hue_range(_data["size"], size_portion)

			# snorm is the normalizer for size
			size_norm = Normalize(vmin=size_norm[0], vmax=size_norm[1])

		# discard s from _scatter_kws and use size in sns.scatterplot
		s = _scatter_kws.pop("s")
		if sizes is None:
			sizes = (min(s, 1), s)

	sns.scatterplot(
		x="x",
		y="y",
		data=_data,
		ax=ax,
		hue="hue",
		size=size,
		sizes=sizes,
		size_norm=size_norm,
		**_scatter_kws,
	)

	# deal with text annotation
	code2label=None
	if text_anno is not None:
		# data
		if isinstance(text_anno, str):
			_data["text_anno"] = _take_data_series(data, text_anno)
		else:
			_data["text_anno"] = text_anno.copy()
		if str(_data["text_anno"].dtype) == "category":
			_data["text_anno"] = _data["text_anno"].cat.remove_unused_categories()

		# text kws
		text_kws = {} if text_kws is None else text_kws
		default_text_kws = dict(
			color='white',  # color for the text, could be a dict, keys are text to be annotated
			fontweight="bold", #fontsize=labelsize,
			bbox=dict(facecolor=palette_dict, # if None, use default color
				boxstyle='round', #ellipse, round
				edgecolor='white', fill=True, linewidth=linewidth, alpha=alpha))
		# coding & id_marker
		text_anno='text_anno'
		if not coding is None and coding!=False:
			if coding == True:
				_data['code'] = _data['hue'].cat.codes #int
			else:
				assert isinstance(coding,str)
				_data["code"] = _take_data_series(data, coding)
				_data=_data.loc[_data['code'].notna()]
				_data["code"]=_data["code"].astype(int)
			_data["code"] = _data["code"].astype("category").cat.remove_unused_categories()
			text_anno='code'
			_data['color'] = _data['hue'].map(palette_dict)
			code2label=_data.loc[:,['code','hue']].drop_duplicates().set_index('code').hue.to_dict()
			_data['code']=_data['code'].astype(str)
			code_colors=_data.loc[:,['code','color']].drop_duplicates().set_index('code').color.to_dict()
			default_text_kws['bbox']['facecolor']=code_colors # background colors for text annotation
			default_text_kws['bbox']['boxstyle'] = 'circle'
		for k in default_text_kws:
			if k !='bbox':
				text_kws.setdefault(k, default_text_kws[k])
			else:
				if 'bbox' not in text_kws:
					text_kws['bbox']={}
				for k1 in default_text_kws['bbox']:
					text_kws['bbox'].setdefault(k1, default_text_kws['bbox'][k1])

		_text_anno_scatter(
			data=_data[["x", "y", text_anno]],
			ax=ax,
			x="x",
			y="y",
			dodge_text=dodge_text,
			dodge_kws=dodge_kws,
			text_transform=text_transform,
			anno_col=text_anno,
			text_kws=text_kws,
			luminance=luminance,
		)

	# deal with outline
	if not outline is None:
		if isinstance(outline, str):
			_data["outline"] = _take_data_series(data, outline)
		else:
			_data["outline"] = outline.copy()
		_outline_kws = {
			"linewidth": linewidth,
			"palette": None,
			"c": "lightgray",
			"single_contour_pad": outline_pad,
		}
		if outline_kws is not None:
			_outline_kws.update(outline_kws)
		density_contour(ax=ax, data=_data, x="x", y="y", groupby="outline", **_outline_kws)

	# clean axis
	if axis_format == "tiny":
		_make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
	elif (axis_format == "empty") or (axis_format is None):
		despine(ax=ax, left=True, bottom=True)
		ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
	else:
		pass

	# deal with legend
	if show_legend and (hue is not None):
		n_hue = len(palette_dict)
		ncol=1 if n_hue <= 40 else 2 if n_hue <= 100 else 3
		if legend_kws is None:
			legend_kws = {}
		default_lgd_kws = dict(
			ncol=ncol,fontsize=labelsize,
			bbox_to_anchor=(1, 1),loc="upper left",
			# borderpad=0.4, # pad between marker (text) and border
			# labelspacing=0.2, #The vertical space between the legend entries, in font-size units
			# handleheight=0.5, #The height of the legend handles, in font-size units.
			# handletextpad=0.2, # The pad between the legend handle (marker) and text, in font-size units.
			# borderaxespad=0.3, # The pad between the Axes and legend border, in font-size units
			# columnspacing=0.2, #The spacing between columns, in font-size units
			markersize=labelsize #legend_kws["fontsize"],
		)
		for k in default_lgd_kws:
			legend_kws.setdefault(k, default_lgd_kws[k])

		exist_hues = _data["hue"].unique()
		color_dict={hue_name: color for hue_name, color in palette_dict.items() if hue_name in exist_hues}
		
		if not code2label is None and id_marker:
			boxstyle='Circle' if not rectangle_marker else 'Round'
			plot_text_legend(color_dict, code2label, ax, title=hue, 
					color_text=legend_color_text, boxstyle=boxstyle,marker_pad=marker_pad,
					legend_kws=legend_kws,marker_fontsize=marker_fontsize,
					alpha=alpha,luminance=luminance)
		else:
			if rectangle_marker:
				## plot Patch legend (rectangle marker)
				plot_color_dict_legend(
					D=color_dict, ax=ax, title=hue, color_text=legend_color_text, 
					kws=legend_kws,luminance=luminance
				)
			else:
				# plot marker legend (for example, circle marker)
				plot_marker_legend(
					color_dict=color_dict, ax=ax, title=hue, color_text=legend_color_text, 
					marker='o',kws=legend_kws,luminance=luminance
				)

	if zoomxy is not None:
		zoom_ax(ax, zoomxy)

	return _data

def get_cmap(cmap):
	try:
		return plt.colormaps.get(cmap)  # matplotlib >= 3.5.1?
	except:
		return plt.get_cmap(cmap)  # matplotlib <=3.4.3?
	
def continuous_scatter(
    data,
    ax=None,
    coord_base="umap",
    x=None,
    y=None,
    scatter_kws=None,
    hue=None,
    hue_norm=None,
    hue_portion=0.95,
    color=None,
    cmap="viridis",
    colorbar=True,
    size=None,
    size_norm=None,
    size_portion=0.95,
    sizes=None,
    sizebar=True,
    text_anno=None,
    dodge_text=False,
    dodge_kws=None,
    text_kws=None,luminance=0.48,
    text_transform=None,
    axis_format="tiny",
    max_points=50000,
    s="auto",
    labelsize=6,
	ticklabel_size=4,
    linewidth=0.5,
    zoomxy=1.05,
    outline=None,
    outline_kws=None,
    outline_pad=2,
    return_fig=False,
    rasterized="auto",
    cbar_kws=None,cbar_width=3,
):
	"""
	Plot scatter on given adata.

	Parameters
	----------
	data : _type_
		_description_
	ax : _type_, optional
		_description_, by default None
	coord_base : str, optional
		_description_, by default "umap"
	x : _type_, optional
		_description_, by default None
	y : _type_, optional
		_description_, by default None
	scatter_kws : _type_, optional
		_description_, by default None
	hue : _type_, optional
		_description_, by default None
	hue_norm : _type_, optional
		_description_, by default None
	hue_portion : float, optional
		_description_, by default 0.95
	color : _type_, optional
		_description_, by default None
	cmap : str, optional
		_description_, by default "viridis"
	colorbar : bool, optional
		_description_, by default True
	size : _type_, optional
		_description_, by default None
	size_norm : _type_, optional
		_description_, by default None
	size_portion : float, optional
		_description_, by default 0.95
	sizes : _type_, optional
		_description_, by default None
	sizebar : bool, optional
		_description_, by default True
	text_anno : _type_, optional
		_description_, by default None
	dodge_text : bool, optional
		_description_, by default False
	dodge_kws : _type_, optional
		_description_, by default None
	text_kws : _type_, optional
		_description_, by default None
	luminance : float, optional
		_description_, by default 0.48
	text_transform : _type_, optional
		_description_, by default None
	axis_format : str, optional
		_description_, by default "tiny"
	max_points : int, optional
		_description_, by default 50000
	s : str, optional
		_description_, by default "auto"
	labelsize : int, optional
		_description_, by default 6
	ticklabel_size : int, optional
		_description_, by default 4
	linewidth : float, optional
		_description_, by default 0.5
	zoomxy : float, optional
		_description_, by default 1.05
	outline : _type_, optional
		_description_, by default None
	outline_kws : _type_, optional
		_description_, by default None
	outline_pad : int, optional
		_description_, by default 2
	return_fig : bool, optional
		_description_, by default False
	rasterized : str, optional
		_description_, by default "auto"
	cbar_kws : _type_, optional
		_description_, by default None
	cbar_width : int, optional
		width of colorbar, by default 3 mm

	Returns
	-------
	_type_
		_description_

	Raises
	------
	ValueError
		_description_
	TypeError
		_description_
	"""
	import seaborn as sns
	import copy
	from matplotlib.cm import ScalarMappable
	# init figure if not provided
	if ax is None:
		fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
	else:
		fig = None

	# add coords
	_data, x, y = _extract_coords(data, coord_base, x, y)
	# _data has 2 cols: "x" and "y"

	# down sample plot data if needed.
	if max_points is not None:
		if _data.shape[0] > max_points:
			_data = _density_based_sample(_data, seed=1, size=max_points, coords=["x", "y"])
	n_dots = _data.shape[0]

	# determine rasterized
	if rasterized == "auto":
		if n_dots > 200:
			rasterized = True
		else:
			rasterized = False

	# auto size if user didn't provide one
	if s == "auto":
		s = _auto_size(ax, n_dots)

	# default scatter options
	_scatter_kws = {"linewidth": 0, "s": s, "legend": None, "rasterized": rasterized}
	if color is not None:
		if hue is not None:
			raise ValueError("Only one of color and hue can be provided")
		_scatter_kws["color"] = color
	if scatter_kws is not None:
		_scatter_kws.update(scatter_kws)

	# deal with color
	if hue is not None:
		if isinstance(hue, str):
			_data["hue"] = _take_data_series(data, hue).astype(float)
			colorbar_label = hue
		else:
			_data["hue"] = hue.astype(float)
			colorbar_label = hue.name

		if hue_norm is None:
			# get the smallest range that include "hue_portion" of data
			# hue_norm = tight_hue_range(_data["hue"], hue_portion)
			hue_norm=(_data["hue"].quantile(1-hue_portion),_data["hue"].quantile(hue_portion))
		# cnorm is the normalizer for color
		cnorm = Normalize(vmin=hue_norm[0], vmax=hue_norm[1])
		if isinstance(cmap, str):
			# from here, cmap become colormap object
			cmap = copy.copy(get_cmap(cmap))
			cmap.set_bad(color=(0.5, 0.5, 0.5, 0.5))
		else:
			if not isinstance(cmap, ScalarMappable):
				raise TypeError(f"cmap can only be str or ScalarMappable, got {type(cmap)}")
	else:
		hue_norm = None
		cnorm = None
		colorbar_label = ""

	# deal with size
	if size is not None:
		if isinstance(size, str):
			_data["size"] = _take_data_series(data, size).astype(float)
		else:
			_data["size"] = size.astype(float)
		size = "size"

		if size_norm is None:
			# get the smallest range that include "size_portion" of data
			size_norm = tight_hue_range(_data["size"], size_portion)

			# snorm is the normalizer for size
			size_norm = Normalize(vmin=size_norm[0], vmax=size_norm[1])

		# replace s with sizes
		s = _scatter_kws.pop("s")
		if sizes is None:
			sizes = (min(s, 1), s)
	else:
		size_norm = None
		sizes = None

	sns.scatterplot(
		x="x",
		y="y",
		data=_data,
		hue="hue",
		palette=cmap,
		hue_norm=cnorm,
		size=size,
		sizes=sizes,
		size_norm=size_norm,
		ax=ax,
		**_scatter_kws,
	)

	if text_anno is not None:
		if isinstance(text_anno, str):
			_data["text_anno"] = _take_data_series(data, text_anno)
		else:
			_data["text_anno"] = text_anno
		if str(_data["text_anno"].dtype) == "category":
			_data["text_anno"] = _data["text_anno"].cat.remove_unused_categories()

		_text_anno_scatter(
			data=_data[["x", "y", "text_anno"]],
			ax=ax,
			x="x",
			y="y",
			dodge_text=dodge_text,
			dodge_kws=dodge_kws,
			text_transform=text_transform,
			anno_col="text_anno",
			text_kws=text_kws,
			luminance=luminance,
		)

	# deal with outline
	if outline:
		if isinstance(outline, str):
			_data["outline"] = _take_data_series(data, outline)
		else:
			_data["outline"] = outline
		_outline_kws = {
			"linewidth": linewidth,
			"palette": None,
			"c": "lightgray",
			"single_contour_pad": outline_pad,
		}
		if outline_kws is not None:
			_outline_kws.update(outline_kws)
		density_contour(ax=ax, data=_data, x="x", y="y", groupby="outline", **_outline_kws)

	# clean axis
	if axis_format == "tiny":
		_make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
	elif (axis_format == "empty") or (axis_format is None):
		despine(ax=ax, left=True, bottom=True)
		ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
	else:
		pass

	return_axes = [ax]

	# make color bar
	if colorbar and (hue is not None):
		# small ax for colorbar
		# default_cbar_kws=dict(loc="upper left", borderpad=0,width="3%", height="20%") #bbox_to_anchor=(1,1)
		if cbar_kws is None:
			cbar_kws={}
		# for k in default_cbar_kws:
		#     if k not in cbar_kws:
		#         cbar_kws[k]=default_cbar_kws[k]

		mm2inch = 1 / 25.4
		space=0
		legend_width = (
			cbar_width * mm2inch * ax.figure.dpi / ax.figure.get_window_extent().width
		)  # mm to px to fraction
		pad = (space + ax.yaxis.labelpad * 1.2 * ax.figure.dpi / 72) / ax.figure.get_window_extent().width
		# labelpad unit is points
		left = ax.get_position().x1 + pad
		ax_legend = ax.figure.add_axes(
			[left, ax.get_position().height * 0.8, legend_width, ax.get_position().height * 0.2]
		)  # left, bottom, width, height
		# print("test:",hue_norm)
		# cbar_kws.setdefault('vmin',hue_norm[0])
		# cbar_kws.setdefault('vmax',hue_norm[1])
		cbar_kws['vmin']=hue_norm[0]
		cbar_kws['vmax']=hue_norm[1]
		cbar = plot_cmap_legend(
			ax=ax,
			cax=ax_legend,
			cmap=cmap,
			label=hue,
			kws=cbar_kws.copy(),labelsize=labelsize, 
			linewidth=linewidth,ticklabel_size=ticklabel_size,
			)
		return_axes.append([ax_legend,cbar])

	# make size bar
	if sizebar and (size is not None):
		# TODO plot dot size bar
		pass

	if zoomxy is not None:
		zoom_ax(ax, zoomxy)

	if return_fig:
		return (fig, tuple(return_axes)), _data
	else:
		return

def plot_cluster(
	adata_path,ax=None,coord_base='tsne',cluster_col='MajorType',
	palette_path=None,coding=True,id_marker=True,
	output=None,
	show=True,figsize=(4, 3.5),sheet_name=None,
	ncol=None,fontsize=5,legend_fontsize=5,
	legend_kws=None,legend_title_fontsize=5,
	marker_fontsize=4,marker_pad=0.1,
	linewidth=0.5,axis_format='tiny',alpha=0.7,
	text_kws=None,**kwargs):
	"""
	Plot cluster.

	Parameters
	----------
	adata_path :
	ax :
	coord_base :
	cluster_col :
	palette_path :
	coding :
	output :
	show :
	figsize :
	sheet_name :
	ncol :
	fontsize :
	legend_fontsize : int
		legend fontsize, default 5
	legend_kws: dict
		kwargs passed to ax.legend
	legend_title_fontsize: int
		legend title fontsize, default 5
	marker_fontsize: int
		Marker fontsize, default 3
		if id_marker is True, and coding is True. legend marker will be a circle (or rectangle) with code
	linewidth : float
		Line width of the legend marker (circle or rectangle), default 0.5
	kwargs : dict
		set text_anno=None to plot clustering without text annotations,
		coding=True to plot clustering without code annotations,
		set show_legend=False to remove the legend

	Returns
	-------

	"""
	from pandas.api.types import is_categorical_dtype
	if coord_base.startswith("X_"):
		coord_base=coord_base.replace('X_','')
	if sheet_name is None:
		sheet_name=cluster_col
	if isinstance(adata_path,str):
		adata=anndata.read_h5ad(adata_path,backed='r')
	else:
		adata=adata_path
	if not is_categorical_dtype(adata.obs[cluster_col]):
		adata.obs[cluster_col] = adata.obs[cluster_col].astype('category')
	if not palette_path is None:
		if isinstance(palette_path,str):
			colors=pd.read_excel(os.path.expanduser(palette_path),sheet_name=sheet_name,index_col=0).Hex.to_dict()
			keys=list(colors.keys())
			existed_vals=adata.obs[cluster_col].unique().tolist()
			for k in existed_vals:
				if k not in keys:
					colors[k]='gray'
			for k in keys:
				if k not in existed_vals:
					del colors[k]
		else:
			colors=palette_path
		adata.uns[cluster_col + '_colors'] = [colors.get(k, 'grey') for k in adata.obs[cluster_col].cat.categories.tolist()]
	else:
		if f'{cluster_col}_colors' not in adata.uns:
			sc.pl.embedding(adata,basis=f"X_{coord_base}",color=[cluster_col],show=False)
		colors={cluster:color for cluster,color in zip(adata.obs[cluster_col].cat.categories.tolist(),adata.uns[f'{cluster_col}_colors'])}

	hue=cluster_col
	text_anno = cluster_col
	text_kws = {} if text_kws is None else text_kws
	text_kws.setdefault("fontsize", fontsize)
	kwargs.setdefault("hue",hue)
	kwargs.setdefault("text_anno", text_anno)
	kwargs.setdefault("text_kws", text_kws)
	kwargs.setdefault("luminance", 0.65)
	kwargs.setdefault("dodge_text", False)
	kwargs.setdefault("axis_format", axis_format)
	kwargs.setdefault("show_legend", True)
	kwargs.setdefault("marker_fontsize", marker_fontsize)
	kwargs.setdefault("marker_pad", marker_pad)
	kwargs.setdefault("linewidth", linewidth)
	kwargs.setdefault("alpha", alpha)
	kwargs["coding"]=coding
	kwargs["id_marker"]=id_marker
	legend_kws={} if legend_kws is None else legend_kws
	default_lgd_kws=dict(
		fontsize=legend_fontsize,
		title=cluster_col,title_fontsize=legend_title_fontsize)
	if not ncol is None:
		default_lgd_kws['ncol']=ncol
	for k in default_lgd_kws:
		legend_kws.setdefault(k, default_lgd_kws[k])
	kwargs.setdefault("dodge_kws", {
			"arrowprops": {
				"arrowstyle": "->",
				"fc": 'grey',
				"ec": "none",
				"connectionstyle": "angle,angleA=-90,angleB=180,rad=5",
			},
			'autoalign': 'xy'})
	if ax is None:
		fig, ax = plt.subplots(figsize=figsize, dpi=300)
	p = categorical_scatter(
		data=adata[adata.obs[cluster_col].notna(),],
		ax=ax,
		coord_base=coord_base,
		palette=colors,legend_kws=legend_kws,
		**kwargs)

	if not output is None:
		plt.savefig(os.path.expanduser(output),transparent=True,bbox_inches='tight',dpi=300)
	if show:
		plt.show()

def plot_gene(
	adata_path="~/Projects/BG/adata/BG.gene-CGN.h5ad",ann_path=None,
	group_col=None,gene='CADM1',query_str=None,title=None,
	palette_path=None,hue_norm=None,
	cbar_kws=dict(extendfrac=0.1),axis_format="tiny",scatter_kws={},
	obsm=None,coord_base='umap',normalize_per_cell=True,
	stripplot=False,hypo_score=False,ylim=None,
	clip_norm_value=10,min_cells=3,cmap='parula',figdir="figures"):
	"""
	Plot gene expression in a given region or group on embedding of adata.

	Parameters
	----------
	adata_path : str, optional
		_description_, by default "~/Projects/BG/adata/BG.gene-CGN.h5ad"
	group_col : str, optional
		_description_, for example: 'Region', by default None
	gene : str, optional
		_description_, by default 'CADM1'
	query_str : _type_, optional
		_description_, by default None
	title : _type_, optional
		_description_, by default None
	palette_path : str, optional
		_description_, by default "~/Projects/BG/metadata/HMBA_color_palette.xlsx"
	obsm : str, optional
		_description_, by default "~/Projects/BG/clustering/100kb/annotated.adata.h5ad"
	coord_base : str, optional
		_description_, by default 'umap'
	normalize_per_cell : bool, optional
		_description_, by default True
	stripplot : bool, optional
		_description_, by default False
	clip_norm_value : int, optional
		_description_, by default 10
	min_cells : int, optional
		_description_, by default 3
	cmap: str, optional
		_description_, by default 'parula'
	figdir : str, optional
		_description_, by default "figures"
	"""

	# sc.set_figure_params(dpi=100,dpi_save=300,frameon=False)
	if title is None:
		if not query_str is None:
			title=query_str
		else:
			title=group_col if not group_col is None else gene
	if not os.path.exists(figdir):
		os.makedirs(figdir, exist_ok=True)
	raw_adata = anndata.read_h5ad(os.path.expanduser(adata_path), backed='r')
	adata = raw_adata[:, gene].to_memory()
	raw_adata.file.close() # close the file to save memory
	if normalize_per_cell:
		adata = normalize_mc_by_cell(
			use_adata=adata, normalize_per_cell=normalize_per_cell,
			clip_norm_value=clip_norm_value,hypo_score=hypo_score)
	is_open=False
	if not obsm is None:
		if isinstance(obsm, str):
			obsm = anndata.read_h5ad(os.path.expanduser(obsm),backed='r')
			is_open=True
		assert isinstance(obsm, anndata.AnnData), "obsm should be an anndata object or a path to an h5ad file."
		keep_cells = list(set(adata.obs_names.tolist()) & set(obsm.obs_names.tolist()))
		adata = adata[keep_cells, :]
		adata.obsm = obsm[keep_cells].obsm
		cur_cols = adata.obs.columns.tolist()
		for col in obsm.obs.columns.tolist():
			if col not in cur_cols:
				adata.obs[col] = obsm.obs.loc[adata.obs_names, col].tolist()
	if is_open:
		obsm.file.close()
	if not ann_path is None:
		if isinstance(ann_path,str):
			metadata=pd.read_csv(os.path.expanduser(ann_path),
				sep='\t',index_col=0)
		else:
			metadata=ann_path.copy()
	else:
		metadata=adata.obs.copy()
	if not query_str is None:
		metadata = metadata.query(query_str)
	overlapped_cells=list(set(adata.obs_names.tolist()) & set(metadata.index.tolist()))
	metadata=metadata.loc[overlapped_cells]
	adata=adata[overlapped_cells,:] # type: ignore
	adata.obs=metadata.loc[adata.obs_names.tolist()]
	print(adata.shape)
	# read color palette
	if not group_col is None and not palette_path is None:
		if os.path.exists(os.path.expanduser(palette_path)):
			palette_path = os.path.abspath(os.path.expanduser(palette_path))
			D = pd.read_excel(palette_path,
							  sheet_name=None, index_col=0)
			color_palette = D[group_col].Hex.to_dict()
		else:
			color_palette = adata.obs.reset_index().loc[:, [group_col, \
				palette_path]].drop_duplicates().dropna().set_index(group_col)[
				palette_path].to_dict()
	else:
		color_palette = None
	# plot gene on given cordinate base
	fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
	output=os.path.join(figdir, f"{title}.{gene}.{coord_base}.pdf")
	sc.pl.embedding(adata, basis=coord_base,
					wspace=0.1, color=[gene],use_raw=False,
					ncols=2, vmin='p5', vmax='p95', frameon=False,
					show=False,cmap=cmap,ax=ax)
	colorbar = fig.axes[-1]
	cur_pos=colorbar.get_position()
	colorbar.set_position([cur_pos.x0,(1-cur_pos.height/2)/2,cur_pos.width, cur_pos.height / 2])
	fig.savefig(output, transparent=True,bbox_inches='tight',dpi=300)

	adata.obs[gene]=adata.to_df().loc[adata.obs_names.tolist(), gene].tolist()
	# print(hue_norm)
	fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
	output=os.path.join(figdir, f"{title}.{gene}.{coord_base}1.pdf")
	continuous_scatter(
		data=adata,
		ax=ax,cmap=cmap,
		hue_norm=hue_norm,
		cbar_kws=cbar_kws,
		hue=gene,axis_format=axis_format,
		text_anno=None,
		coord_base=coord_base,**scatter_kws)
	fig.savefig(output, transparent=True,bbox_inches='tight',dpi=300)
	
	if not group_col is None:
		output=os.path.join(figdir, f"{title}.{group_col}.{coord_base}.pdf")
		if not os.path.exists(output): # plot embedding colored by group_col
			if not color_palette is None:
				use_cells = adata.obs.loc[adata.obs[group_col].isin(list(color_palette.keys()))].index.tolist()
			else:
				use_cells = adata.obs.index.tolist()
			plot_cluster(adata_path=adata,coord_base=coord_base,
				cluster_col=group_col,
				coding=False,palette_path=palette_path,ncol=1,
				output=output,text_anno=None)

		# boxplot
		data = adata.to_df()
		data[group_col] = adata.obs.loc[data.index.tolist(), group_col].tolist()
		vc = data[group_col].value_counts()
		N=vc.shape[0]
		if not color_palette is None:
			keep_groups = list(set(list(color_palette.keys())) & set(vc[vc >= min_cells].index.tolist()))
			data = data.loc[data[group_col].isin(keep_groups)]
		vc = vc.to_dict()
		order = data.groupby(group_col)[gene].median().sort_values().index.tolist()
		width = max(5, N*0.5)
		plt.figure(figsize=(width, 3.5))
		if stripplot:
			ax = sns.stripplot(data=data, jitter=0.4,
							edgecolor='white', x=group_col, y=gene, palette=color_palette, \
							order=order, size=0.5)
		else:
			ax = None
		# ax = sns.boxplot(data=data, x=group_col, y=gene, palette=color_palette, ax=ax,  # hue=group_col,
		# 				fliersize=0.5, notch=False, showfliers=False, saturation=0.6, order=order)
		# boxplot are incorrect for some cases when there are many 0, median and lower quartile are often at zero; use violinplot in stead.
		ax = sns.violinplot(data=data, x=group_col, y=gene, palette=color_palette, ax=ax,  # hue=group_col,
                saturation=0.6, order=order,density_norm='width',cut=0,bw_adjust=0.5)
		# ax=sns.swarmplot(data=data,palette=color_palette,\
		#                   edgecolor='white',x=group_col,y=gene,\
		#                   order=order)
		if not ylim is None:
			ax.set_ybound(ylim)
		ax.set_xticklabels([f"{label} ({vc[label]})" for label in order])
		title=title.replace(' ','.')
		ax.set_title(title)
		ax.xaxis.label.set_visible(False)
		plt.setp(ax.xaxis.get_majorticklabels(), rotation=-45, ha='left')
		plt.savefig(os.path.join(figdir, f"{title}.{gene}.{group_col}.boxplot.pdf"), transparent=True,bbox_inches='tight',dpi=300)
	return adata

def plot_genes(
	adata_path="/home/x-wding2/Projects/BICAN/adata/HMBA_v2/HMBA.Group.downsample_1500.h5ad",
	query_str=None,
	ann_path=None, #"~/Projects/BG/clustering/100kb/annotations.tsv",
	group_col='Subclass',
	parent_col=None,
	modality='RNA', # mc or RNA
	use_raw=True, # True for RNA
	expression_cutoff='p5', # for RNA, could be int, median, mean of p5, p95 and so on
	genes=None,
	cell_type_order=None,
	gene_order=None,
	row_cluster=False,
	col_cluster=False,
	cmap='Greens_r',
	group_legend=False,
	parent_legend=False,
	title=None,
	palette_path=None,#"/home/x-wding2/Projects/BICAN/adata/HMBA_v2/HMBA_color_palette.xlsx"
	legend_kws=dict(extendfrac=0.1,extend='both',label='Mean mCG'),
	normalize_per_cell=True,
	clip_norm_value=10,
	hypo_score=False,
	figsize=(10, 3.5),
	marker='o',
	plot_kws={},transpose=False,
	outname="test.pdf"):
	"""
	_summary_

	Parameters
	----------
	couldbeint : _type_
		_description_
	median : _type_
		_description_
	meanofp5 : _type_
		_description_
	adata_path : str, optional
		_description_, by default "/home/x-wding2/Projects/BICAN/adata/HMBA_v2/HMBA.Group.downsample_1500.h5ad"
	query_str : _type_, optional
		_description_, by default None
	ann_path : _type_, optional
		_description_, by default None
	group_col : str, optional
		_description_, by default 'Subclass'
	parent_col : str, optional
		_description_, by default 'Neighborhood'
	modality : str, optional
		_description_, by default 'RNA'
	p95andsoongenes : _type_, optional
		_description_, by default None
	cell_type_order : _type_, optional
		_description_, by default None
	gene_order : _type_, optional
		_description_, by default None
	row_cluster : bool, optional
		_description_, by default False
	col_cluster : bool, optional
		_description_, by default False
	cmap : str, optional
		_description_, by default 'Greens_r'
	group_legend : bool, optional
		_description_, by default False
	parent_legend : bool, optional
		_description_, by default False
	title : str, optional
		_description_, by default 'test'
	palette_path : _type_, optional
		_description_, by default None
	obsm : _type_, optional
		_description_, by default None
	normalize_per_cell : bool, optional
		_description_, by default True
	clip_norm_value : int, optional
		_description_, by default 10
	hypo_score : bool, optional
		_description_, by default False
	figsize : tuple, optional
		_description_, by default (10, 3.5)
	cmap : str, optional
		_description_, by default 'Greens_r'
	marker : str, optional
		_description_, by default 'o'
	plot_kws : dict, optional
		_description_, by default {}
	outname : str, optional
		_description_, by default "test.pdf"
	"""
	assert not genes is None, "Please provide genes to plot."
	# adata could be single cell level or pseudobulk level (adata.layers['frac'] should be existed)
	raw_adata = anndata.read_h5ad(os.path.expanduser(adata_path), backed='r')

	all_vars=set(raw_adata.var_names.tolist())
	keep_genes=list(set(all_vars) & set(genes)) # keep_genes=[g for g in all_vars if g in genes]
	error_genes=[g for g in genes if g not in keep_genes]
	if len(error_genes)>0:
		print(f"genes not found in adata: {error_genes}")
	adata = raw_adata[:, keep_genes].to_memory() # type: ignore
	if use_raw and not adata.raw is None:
		adata_raw=adata.raw[:,adata.var_names.tolist()].to_adata()
		adata.X=adata_raw[adata.obs_names.tolist(),adata.var_names.tolist()].X.copy() # type: ignore
		del adata_raw
	raw_adata.file.close() # close the file to save memory

	if not ann_path is None:
		if isinstance(ann_path,str):
			metadata=pd.read_csv(os.path.expanduser(ann_path),
				sep='\t',index_col=0)
		else:
			metadata=ann_path.copy()
	else:
		metadata=adata.obs.copy()
	if not query_str is None:
		metadata = metadata.query(query_str)
	overlapped_cells=list(set(adata.obs_names.tolist()) & set(metadata.index.tolist()))
	metadata=metadata.loc[overlapped_cells]
	adata=adata[overlapped_cells,:] # type: ignore
	if isinstance(group_col,list):
		group_col1="+".join(group_col)
		metadata[group_col1]=metadata.loc[:,group_col].apply(lambda x:'+'.join(x.astype(str).tolist()),axis=1)
		group_col=group_col1
	adata.obs[group_col]=metadata.loc[adata.obs_names.tolist(),group_col].tolist()
	if title is None:
		if not query_str is None:
			title=query_str
		else:
			title=group_col if not group_col is None else '-'.join(genes)
	if not parent_col is None and parent_col not in adata.obs.columns.tolist():
		adata.obs[parent_col]=metadata.loc[adata.obs_names.tolist(),parent_col].tolist()
			
	if modality!='RNA' and normalize_per_cell:
		adata = normalize_mc_by_cell(
			use_adata=adata, normalize_per_cell=normalize_per_cell,
			clip_norm_value=clip_norm_value,hypo_score=hypo_score)
	print(adata.shape)

	# read color palette
	color_palette={}
	if not palette_path is None:
		if os.path.exists(os.path.expanduser(palette_path)):
			palette_path = os.path.abspath(os.path.expanduser(palette_path))
			D = pd.read_excel(palette_path,
							sheet_name=None, index_col=0)
			if group_col in D:
				color_palette[group_col] = D[group_col].Hex.to_dict()
			else:
				assert '+' in group_col, f"{group_col} not found in the palette file."
				for group in group_col.split('+'):
					assert group in D, f"{group} not found in the palette file."
					color_palette[group] = D[group].Hex.to_dict()
			if not parent_col is None:
				color_palette[parent_col] = D[parent_col].Hex.to_dict()
		else:
			color_palette[group_col] = adata.obs.reset_index().loc[:, [group_col, \
				palette_path]].drop_duplicates().dropna().set_index(group_col)[
				palette_path].to_dict()
			color_palette[parent_col] = adata.obs.reset_index().loc[:, [parent_col, \
				palette_path]].drop_duplicates().dropna().set_index(parent_col)[
				palette_path].to_dict()
	else:
		color_palette = None

	data=adata.to_df() # rows are cells or cell types, columns are genes
	if modality=='RNA' and isinstance(expression_cutoff,str):
		if expression_cutoff=='median':
			cutoff=data.stack().median()
		elif expression_cutoff=='mean':
			cutoff=data.stack().mean()
		else: # quantile, such as p5,p95
			f=float(expression_cutoff.replace('p',''))
			cutoff=data.stack().quantile(f/100)
		expression_cutoff=cutoff
			
	data[group_col]=adata.obs.loc[data.index.tolist(),group_col].tolist()
	if not parent_col is None and parent_col in adata.obs.columns.tolist():
		group2parent=adata.obs.loc[:,[group_col,parent_col]].drop_duplicates().set_index(group_col)[parent_col].to_dict()
	plot_data=data.groupby(group_col).mean().stack().reset_index()
	plot_data.columns=[group_col,'Gene','Mean']
	if 'frac' in adata.layers:
		D=adata.to_df(layer='frac').stack().to_dict()
	else:
		if modality!='RNA': # methylation, cutoff = 1
			assert normalize_per_cell==True,"Normalized methylation fraction is required"
			hypo_frac=data.groupby(group_col).agg(lambda x:x[x< 1].shape[0] / x.shape[0]) # fraction of cells showing hypomethylation for the corresponding genes
			D=hypo_frac.stack().to_dict()
		else: # for RNA
			print(f"Using expression cutoff: {expression_cutoff}")
			frac=data.groupby(group_col).agg(lambda x:x[x>expression_cutoff].shape[0] / x.shape[0]) # raw count > expression_cutoff means the gene is expressed
			D=frac.stack().to_dict()
	plot_data['frac']=plot_data.loc[:,[group_col,'Gene']].apply(lambda x:tuple(x.tolist()),axis=1).map(D)
	# plot_data

	df_cols=pd.DataFrame(list(sorted(adata.obs[group_col].unique().tolist())),columns=[group_col])
	if not parent_col is None:
		df_cols[parent_col]=df_cols[group_col].map(group2parent)
		df_cols.sort_values([parent_col,group_col],inplace=True)
	df_cols.index=df_cols[group_col].tolist()
	if not cell_type_order is None:
		rows=[ct for ct in cell_type_order if ct in df_cols.index.tolist()]
		df_cols=df_cols.loc[rows]
	col_ha_dict={}
	if '+' in group_col:
		individual_groups=group_col.split('+')
		for ig in individual_groups:
			df_cols[ig]=df_cols[group_col].apply(lambda x:x.split('+')[individual_groups.index(ig)])
			group_colors={}
			for k in df_cols[ig].unique().tolist():
				group_colors[k]=color_palette[ig][k]
			col_ha_dict[ig]=anno_simple(df_cols[ig],colors=group_colors,
								add_text=False,legend=group_legend,height=3,label=ig)
	df_cols.dropna(inplace=True)
	# df_cols.head()
	if not parent_col is None:
		parent_colors={}
		axis=1 if not transpose else 0 # 1 for vertical (col annotation), 0 for horizontal
		for k in df_cols[parent_col].unique().tolist():
			parent_colors[k]=color_palette[parent_col][k]
		if '+' not in group_col:
			group_colors={}
			for k in df_cols[group_col].unique().tolist():
				group_colors[k]=color_palette[group_col][k]
			col_ha=HeatmapAnnotation(axis=axis,
				label=anno_label(df_cols[group_col], colors=group_colors,merge=True,
								rotation=45,fontsize=12,arrowprops = dict(visible=False)),
				group=anno_simple(df_cols[group_col],colors=group_colors,
									add_text=False,legend=group_legend,height=3,label=group_col), 
				parent=anno_simple(df_cols[parent_col],colors=parent_colors,
									add_text=False,legend=parent_legend,height=3,label=parent_col), 
			)
		else:
			col_ha_dict[parent_col]=anno_simple(df_cols[parent_col],colors=parent_colors,
									add_text=False,legend=parent_legend,height=3,label=parent_col)
			col_ha = HeatmapAnnotation(**col_ha_dict,axis=axis,
									verbose=0)
		colnames=False
		
	else:
		axis=1 if not transpose else 0 # 1 for vertical (col annotation), 0 for horizontal
		if '+' not in group_col:
			group_colors={}
			for k in df_cols[group_col].unique().tolist():
				group_colors[k]=color_palette[group_col][k]
			col_ha=HeatmapAnnotation(axis=axis,
				group=anno_simple(df_cols[group_col],colors=group_colors,
									add_text=False,legend=group_legend,height=3,label=group_col), 
			)
		else:
			col_ha = HeatmapAnnotation(**col_ha_dict,axis=axis,
									verbose=0)
		colnames=True
	if not transpose:
		top_annotation=col_ha
		left_annotation=None
		x=group_col
		y='Gene'
		x_order=df_cols.index.tolist()
		y_order=gene_order
		show_colnames=colnames
		show_rownames=True
	else:
		top_annotation=None
		left_annotation=col_ha
		y=group_col
		x='Gene'
		y_order=df_cols.index.tolist()
		x_order=gene_order
		show_rownames=colnames
		show_colnames=True

	default_plot_kws=dict(
		marker=marker,grid=None,legend_gap=8,dot_legend_marker=marker,cmap_legend_kws=legend_kws,
		row_cluster=row_cluster,col_cluster=col_cluster,
		row_cluster_method='ward',row_cluster_metric='euclidean',
		col_cluster_method='ward',col_cluster_metric='euclidean',
		col_names_side='top',row_names_side='left',
		show_rownames=show_rownames,show_colnames=show_colnames,row_dendrogram=False,
		# vmin=0,vmax=1.5,
		xticklabels_kws={'labelrotation': 45, 'labelcolor': 'blue','labelsize':10,'top':True},
		yticklabels_kws={'labelcolor': 'blue','labelsize':10,'left':True},
		spines=False,
	)
	for k in default_plot_kws:
		if k not in plot_kws:
			plot_kws[k]=default_plot_kws[k]
			
	plt.figure(figsize=figsize)
	cm1 = DotClustermapPlotter(
		data=plot_data, top_annotation=top_annotation,left_annotation=left_annotation,
		x_order=x_order,y_order=y_order,
		x=x,y=y,value='Mean',c='Mean',s='frac',
		cmap=cmap,verbose=1,**plot_kws,
	)
	for ax in cm1.heatmap_axes.ravel():
		ax.grid(axis='both', which='minor', color='grey', linestyle='--',alpha=0.6,zorder=0)
	if outname is None:
		outname=f"{title}.pdf"
	plt.savefig(os.path.expanduser(outname),transparent=True, bbox_inches='tight',dpi=300)
	plt.show()
	return plot_data,df_cols,cm1