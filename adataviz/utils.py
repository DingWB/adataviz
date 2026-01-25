import pandas as pd
import os, sys
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import matplotlib
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from matplotlib.text import Text
import copy
import warnings
warnings.filterwarnings("ignore")

mm2inch = 1 / 25.4

def df2stdout(df):
	sys.stdout.write('\t'.join([str(i) for i in df.columns.tolist()]) + '\n')
	for i, row in df.iterrows():
		try:
			sys.stdout.write('\t'.join([str(i) for i in row.fillna('').tolist()]) + '\n')
		except:
			sys.stdout.close()

def serialize(x):
	if isinstance(x, pd.DataFrame):
		df2stdout(x)
	elif not x is None:
		print(x)
	else:
		pass

def add_gene_type(input,ref,ref_gene=5,
				  ref_type=6,query_col=6):
	df=pd.read_csv(os.path.abspath(os.path.expanduser(ref)),
				   sep='\t',usecols=[ref_gene,ref_type],header=None)
	df[ref_gene] = df[ref_gene].apply(lambda x: x.capitalize())
	D=df.set_index(ref_gene)[ref_type].to_dict()
	f=open(os.path.abspath(os.path.expanduser(input)),'r')
	for line in f.readlines():
		values=line.rstrip('\n').split('\t')
		gene=values[query_col].capitalize()
		if values[query_col] not in D:
			gene_type=''
		else:
			gene_type=D[values[query_col]]
		values=values+[gene_type]
		try:
			sys.stdout.write('\t'.join(values)+'\n')
		except:
			sys.stdout.close()
			f.close()
			return
	f.close()

def parse_gtf(
	input="gencode.v43.annotation.gtf",
	gene_types=['protein_coding','lncRNA'],
	transcript_types=['protein_coding','lncRNA'],
	gene_distance=[2000,2000],promoter_distance=[1500,2000],
	tabix=None,bgzip=None
):
	if tabix is None:
		tabix=os.path.join(os.path.dirname(sys.executable),"tabix")
	if bgzip is None:
		bgzip=os.path.join(os.path.dirname(sys.executable),"bgzip")
	outdir=os.path.dirname(os.path.abspath(os.path.expanduser(input)))
	gene_up,gene_down=gene_distance
	promoter_up, promoter_down = promoter_distance
	if not os.path.exists(os.path.join(outdir,"parsed_gtf.txt")):
		df=pd.read_csv(os.path.expanduser(input),sep='\t',header=None,
					   comment="#",usecols=[0,2,3,4,6,8],
					   names=['chrom','record_type','beg','end','strand','information'])
		cols=['gene_id','gene_type','gene_name','transcript_id','transcript_type','hgnc_id']
		def parse_info(x):
			x=x.replace('"','')
			D={}
			for item in x.strip().rstrip(';').split(';'):
				try:
					k,v=item.strip().split(' ')
				except:
					print(x,item)
				D[k.strip()]=v.strip()
			return D

		df['info_dict']=df.information.apply(parse_info)
		for col in cols:
			df[col]=df.info_dict.apply(lambda x:x.get(col,''))
		df.hgnc_id=df.hgnc_id.apply(lambda x:x.split(':')[1].strip() if x!='' else x)
		df.drop(['information','info_dict'],axis=1,inplace=True)
		df.to_csv(os.path.join(outdir,"parsed_gtf.txt"),sep='\t',index=False)
	else:
		df=pd.read_csv(os.path.join(outdir,"parsed_gtf.txt"),sep='\t',
					   dtype={'hgnc_id':str})
	# gene body 2k bed
	gdf=df.loc[df.record_type=='gene',['chrom','beg','end','gene_id','gene_name','strand','gene_type']]
	beg = gdf.loc[:, ['beg', 'end', 'strand']].apply(
		lambda x: x.beg - gene_up if x.strand == '+' else x.beg - gene_down, axis=1)
	end = gdf.loc[:, ['beg', 'end', 'strand']].apply(
		lambda x: x.end + gene_down if x.strand == '+' else x.end + gene_up, axis=1)
	gdf.beg = beg.tolist()
	gdf.end = end.tolist()

	if not gene_types is None and len(gene_types) > 0:
		print(gene_types,type(gene_types))
		gdf=gdf.loc[gdf.gene_type.isin(gene_types)]
	print(gdf.shape, gdf.gene_id.nunique(), gdf.iloc[:, :3].drop_duplicates().shape)
	dup = gdf.iloc[:, :3].duplicated()
	if dup.sum() > 0:
		gdf=gdf.groupby(['chrom','beg','end'],as_index=False).agg(
			lambda x:';'.join(set(x))
		)
	print(gdf.shape, gdf.gene_id.nunique(), gdf.iloc[:, :3].drop_duplicates().shape)
	gdf.beg=gdf.beg-1 #bed format, beg is 0-based, end is 1-based
	gdf.gene_type = gdf.gene_type.apply(lambda x: 'protein_coding' if 'protein_coding' in x else x)
	gdf.beg=gdf.beg.apply(lambda x:x if x >=0 else 0)
	outname=os.path.join(outdir,f"gene_{gene_up}_{gene_down}.bed")
	gdf.to_csv(outname,sep='\t',index=False)
	cmd = f"sed '1d' {outname} | sort -k 1,1 -k 2,2n -k 3,3n | {bgzip} >  {outname}.gz"
	print(cmd)
	os.system(cmd)
	cmd = f"{tabix} -0 -b 2 -e 3 {outname}.gz"
	print(cmd)
	os.system(cmd)

	# gene promoter, per gene
	gdf=df.loc[df.record_type=='gene',['chrom','beg','end','gene_id','gene_name','strand','gene_type']]
	if not gene_types is None and len(gene_types) > 0:
		gdf = gdf.loc[gdf.gene_type.isin(gene_types)]
	beg = gdf.loc[:, ['beg', 'end', 'strand']].apply(
		lambda x: x.beg - promoter_up if x.strand == '+' else x.end - promoter_down, axis=1)
	end = gdf.loc[:, ['beg', 'end', 'strand']].apply(
		lambda x: x.beg + promoter_down if x.strand == '+' else x.end + promoter_up, axis=1)
	gdf.beg=beg.tolist()
	gdf.end=end.tolist()
	dup = gdf.iloc[:, :3].duplicated()
	if dup.sum() > 0:
		gdf = gdf.groupby(['chrom', 'beg', 'end'], as_index=False).agg(
			lambda x: ';'.join(set(x))
		)
	gdf.beg = gdf.beg - 1
	gdf.beg = gdf.beg.apply(lambda x: x if x >= 0 else 0)
	gdf.gene_type = gdf.gene_type.apply(lambda x: 'protein_coding' if 'protein_coding' in x else x)
	# write gene promoter to file
	outname = os.path.join(outdir, f"promoter_{promoter_up}_{promoter_down}.bed")
	gdf.to_csv(outname, sep='\t', index=False)
	cmd = f"sed '1d' {outname} | sort -k 1,1 -k 2,2n -k 3,3n | {bgzip} >  {outname}.gz"
	print(cmd)
	os.system(cmd)
	cmd = f"{tabix} -0 -b 2 -e 3 {outname}.gz"
	print(cmd)
	os.system(cmd)

	"""
	# TSS, per transcript
	tdf = df.loc[df.record_type == 'transcript',['chrom', 'beg', 'end', 'transcript_id',
										'gene_name', 'strand', 'transcript_type']]
	if not transcript_types is None and len(transcript_types) > 0:
		tdf=tdf.loc[tdf.transcript_type.isin(transcript_types)]
	beg=tdf.loc[:,['beg', 'end', 'strand']].apply(
		lambda x: x.beg-promoter_up if x.strand == '+' else x.end-promoter_down,axis=1)
	end=tdf.loc[:,['beg', 'end', 'strand']].apply(
		lambda x: x.beg+promoter_down if x.strand == '+' else x.end+promoter_up,axis=1)
	tdf.beg=beg.tolist()
	tdf.end=end.tolist()
	print(tdf.shape, tdf.transcript_id.nunique(), tdf.iloc[:, :3].drop_duplicates().shape)
	dup = tdf.iloc[:, :3].duplicated()
	if dup.sum() > 0:
		tdf = tdf.groupby(['chrom', 'beg', 'end'], as_index=False).agg(
			lambda x: ';'.join(set(x))
		)
	print(tdf.shape, tdf.transcript_id.nunique(), tdf.iloc[:, :3].drop_duplicates().shape)
	tdf.beg = tdf.beg.apply(lambda x: x if x >= 0 else 0)
	# write transcripts promoter to file
	outname = os.path.join(outdir, f"tss_{promoter_up}_{promoter_down}.bed")
	tdf.to_csv(outname, sep='\t', index=False)
	cmd = f"sed '1d' {outname} | sort -k 1,1 -k 2,2n -k 3,3n | {bgzip} >  {outname}.gz"
	print(cmd)
	os.system(cmd)
	cmd = f"{tabix} -0 -b 2 -e 3 {outname}.gz"
	print(cmd)
	os.system(cmd) """

def read_markers(
	marker_path=None,marker_sheet=None,
	max_depth=None,marker_level=None,parent_ct=None,
	species = 'human',marker_cols = ['RNA_markers', 'Methylation_markers']
):
	def markers_to_list(x):
		if not pd.isna(x) and x != '':
			try:
				L = eval(x)
			except:
				L = [g.strip() for g in x.split(',')]
		else:
			L = []
		return L

	def get_unique_genes(x):
		genes=[]
		for g in x:
			if g not in genes:
				genes.append(g)
		return genes

	if marker_path is None:
		return None
	if marker_sheet is None:
		dfs = []
		df = pd.read_excel(marker_path, sheet_name=None)
		for key in df:
			df1 = df[key]
			dfs.append(df1)
		df_marker = pd.concat(dfs)
	else:
		print(marker_sheet)
		df_marker = pd.read_excel(marker_path, sheet_name=marker_sheet)
	for col in marker_cols:
		df_marker[col] = df_marker[col].apply(markers_to_list)
	df_marker['Markers'] = df_marker.loc[:, marker_cols].sum(axis=1).apply(get_unique_genes)
	df_marker = df_marker.loc[df_marker.Markers.apply(lambda x: len(x) > 1)]
	if species.lower() == 'mouse':
		df_marker.Markers = df_marker.Markers.apply(lambda x: [g.capitalize() for g in x])
	else:
		df_marker.Markers = df_marker.Markers.apply(lambda x: [g.upper() for g in x])
	# Filter marker genes
	if max_depth is None:
		if not marker_level is None and marker_level in df_marker.Level.tolist():  # clustering_name_annot: CellClass (L1_annot), MajorType (L2_annot), SubType (L3_annot),
			df_marker = df_marker.loc[df_marker.Level == marker_level]
		if not parent_ct is None and parent_ct in df_marker.Parent.tolist():
			df_marker = df_marker.loc[df_marker.Parent == parent_ct]
	else:
		df_marker = df_marker.loc[df_marker.Depth <= max_depth]
	return df_marker

def summarize_gene_counts(input="gencode.vM23.annotation.gtf",
						  gene_types=['protein_coding','lncRNA']):
	outdir = os.path.dirname(os.path.abspath(os.path.expanduser(input)))
	if not os.path.exists(os.path.join(outdir,"parsed_gtf.txt")):
		df=pd.read_csv(os.path.expanduser(input),sep='\t',header=None,
					   comment="#",usecols=[0,2,3,4,6,8],
					   names=['chrom','record_type','beg','end','strand','information'])
		cols=['gene_id','gene_type','gene_name','transcript_id','transcript_type','hgnc_id']
		def parse_info(x):
			x=x.replace('"','')
			D={}
			for item in x.strip().rstrip(';').split(';'):
				try:
					k,v=item.strip().split(' ')
				except:
					print(x,item)
				D[k.strip()]=v.strip()
			return D

		df['info_dict']=df.information.apply(parse_info)
		for col in cols:
			df[col]=df.info_dict.apply(lambda x:x.get(col,''))
		df.hgnc_id=df.hgnc_id.apply(lambda x:x.split(':')[1].strip() if x!='' else x)
		df.drop(['information','info_dict'],axis=1,inplace=True)
		df.to_csv(os.path.join(outdir,"parsed_gtf.txt"),sep='\t',index=False)
	else:
		df=pd.read_csv(os.path.join(outdir,"parsed_gtf.txt"),sep='\t',
					   dtype={'hgnc_id':str})
	df=df.loc[df.gene_type.isin(gene_types)]
	df=df.loc[df.record_type=='gene']
	vc=df.groupby('chrom').gene_name.nunique()
	vc.sort_values(ascending=False,inplace=True)
	print(vc.to_dict())

def prepare_color_palette(color_dict=None,outpath="palette.xlsx"):
	"""
	Generating a .xlsx file including all color palette.

	Parameters
	----------
	colors : dict
		A dict of dict, keys are categorical terms, values are HEX color code

	Returns
	-------

	"""
	outpath=os.path.expanduser(outpath)
	writer = pd.ExcelWriter(outpath)
	for key in color_dict:
		data = pd.DataFrame.from_dict(color_dict[key], orient='index', columns=['Hex'])
		# data.style.background_gradient(cmap='gray_r')
		# data.style.applymap(lambda x:'color:'+x if x.startswith('#') else 'color: white')
		data.to_excel(writer, sheet_name=key, index=True)
		workbook = writer.book
		worksheet = writer.sheets[key]
		colors = data.Hex.tolist()
		for i in range(data.shape[0]):
			color = colors[i]
			f = workbook.add_format({'bold': True, 'font_color': 'black', 'bg_color': color})
			worksheet.write(i + 1, 1, color, f)
		width = 20
		cell_fmt = workbook.add_format(
			{'bold': False, 'font_color': 'black',
			 # 'bg_color':'green',
			 'align': 'center', 'valign': 'vcenter'})
		# styled = data.style.applymap(lambda val: 'color: %s' % 'red' if val < 0 else 'black').highlight_max()
		worksheet.set_column(0, 1, width, cell_fmt)
	# worksheet.conditional_format(f'A:{last_col}', {'type': 'no_blanks', 'format': cell_fmt})
	writer.close()

def preprocessing_HMBA(
	workdir="~/Projects/BICAN/adata/HMBA_v2",
	adata_path = "~/Projects/BICAN/adata/HMBA_v2/Human_HMBA_basalganglia_AIT_pre-print.h5ad",
	gene_bed="~/Ref/hg38/annotations/gene_2000_2000.bed.gz",
	levels = ['Neighborhood', 'Class', 'Subclass', 'Group'],
	downsample = 1500,topn = 50
):
	import pandas as pd
	import anndata
	import scanpy as sc
	workdir=os.path.expanduser(workdir)
	outdir = os.path.join(workdir,"DEG")
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	rna_cell_type=levels[-1]
	adata = anndata.read_h5ad(os.path.expanduser(adata_path), backed='r')

	# split into Neu and NonN h5ad files
	if not os.path.exists(os.path.join(workdir, "NonN.h5ad")) and not os.path.exists(os.path.join(workdir, "Neu.h5ad")):
		use_cells = adata.obs.query("Neighborhood=='Nonneuron'").index.tolist()
		use_adata = adata[use_cells, :].to_memory()
		use_adata.write_h5ad(os.path.join(workdir, "NonN.h5ad"))

		use_cells = adata.obs.query("Neighborhood!='Nonneuron'").index.tolist()
		use_adata = adata[use_cells, :].to_memory()
		use_adata.write_h5ad(os.path.join(workdir, "Neu.h5ad"))

	# downsample & DEG
	if not os.path.exists(os.path.join(workdir,f"HMBA.{rna_cell_type}.downsample_{downsample}.h5ad")):
		keep_cells = adata.obs.loc[adata.obs[rna_cell_type].notna()].groupby(rna_cell_type).apply(
			lambda x: x.sample(downsample).index.tolist() if x.shape[0] > downsample else x.index.tolist()).sum()
		adata[keep_cells,:].write_h5ad(os.path.join(workdir,f"HMBA.{rna_cell_type}.downsample_{downsample}.h5ad"),
									   compression='gzip')
	adata.file.close()
	rna_adata=anndata.read_h5ad(os.path.join(workdir,f"HMBA.{rna_cell_type}.downsample_{downsample}.h5ad"),backed='r')
	print(rna_adata.obs[rna_cell_type].value_counts())

	# rna_adata.obs.loc[:, levels]
	df_gene = pd.read_csv(gene_bed, sep='\t', header=None)
	df_gene = df_gene.loc[df_gene[6].isin(['protein_coding'])]  # 'lncRNA',

	def get_markers(adata1, level, df_gene, key, outdir,topn,downsample):
		if not os.path.exists(os.path.join(outdir, f"{key}.tsv")):
			use_cells=adata1.obs.groupby(level).apply(
		lambda x: x.sample(downsample).index.tolist() if x.shape[0] > downsample else x.index.tolist()).sum()
			use_adata=adata1[use_cells,:].to_memory()
			vc=use_adata.obs[level].value_counts()
			groups=vc[vc >= 3].index.tolist()
			sc.tl.rank_genes_groups(use_adata, groupby=level, groups=groups,method="wilcoxon", use_raw=False, key_added=key)
			# The Wilcoxon rank-sum test is a non-parametric test that compares the ranks of gene expression values between groups. It works best with raw counts because normalization can alter the rank distribution
			markers = sc.get.rank_genes_groups_df(use_adata, group=groups, key=key)
			markers = markers.loc[~ markers.names.isna()]
			markers = markers.loc[(~ markers.logfoldchanges.isna()) & (markers.scores > 0) & (markers.pvals < 0.01) & (
					markers.logfoldchanges > 1)]
			if markers.shape[0] == 0:
				return {}
			markers.sort_values('logfoldchanges', ascending=False, inplace=True)
			markers.to_csv(os.path.join(outdir, f"{key}.tsv"), sep='\t', index=False)
		else:
			markers=pd.read_csv(os.path.join(outdir, f"{key}.tsv"), sep='\t')
		markers.names = markers.names.apply(lambda x: x.split('.')[0])
		markers = markers.loc[markers.names.isin(df_gene[4].tolist())].groupby('group',
			   include_groups=False).apply(lambda x: x.head(topn).names.tolist()).to_dict()
		return markers

	# DEG
	R = []
	for i in range(len(levels)):
		level = levels[i]
		print(level)
		if i == 0:
			parent = ''
			markers = get_markers(rna_adata, level, df_gene, level,outdir,topn,downsample)
			for k in markers: # k is cell type in different levels
				R.append([level, parent, k, markers[k]]) # append topn marker genes
		else:  # i >= 1
			for clusters, df1 in rna_adata.obs.groupby(levels[:i],observed=True):
				if df1.shape[0]==0:
					continue
				if df1[level].nunique() < 2:
					continue
				print(clusters)
				adata1 = rna_adata[df1.index.tolist(),:].to_memory()
				parent = list(clusters)[-1]  # '|'.join(list(clusters))
				key = parent + '|' + level
				# print(key)
				markers = get_markers(adata1, level, df_gene, key, outdir,topn,downsample)
				if len(markers) == 0:
					continue
				for k in markers:
					R.append([level, parent, k, markers[k]])

	data = pd.DataFrame(R, columns=['Level', 'Parent', 'CellType', 'Markers'])
	data.to_excel(os.path.join(workdir,"HMBA_markers.xlsx"),
				  index=False, sheet_name='HMBA')

	# color palette
	color_dict = {}
	for col in ['Neighborhood', 'Class', 'Subclass', 'Group']:
		D = rna_adata.obs.reset_index().loc[:, [col, f"color_hex_{col.lower()}"]].drop_duplicates().set_index(col)[
			f"color_hex_{col.lower()}"].to_dict()
		color_dict[col] = D
	prepare_color_palette(color_dict=color_dict,
						  outpath=os.path.join(workdir,"HMBA_color_palette.xlsx"))

	adata = anndata.read_h5ad(os.path.expanduser(adata_path), backed='r')
	# add embedding (umap and tsne) to obs
	obs=adata.obs.copy()
	obs.index.name='cell'
	for coord in ["umap",'tsne']:
		if f'X_{coord}' not in adata.obsm:
			continue
		df_coord=pd.DataFrame(adata.obsm[f'X_{coord}'],columns=[f'{coord}_0',f'{coord}_1'],index=adata.obs_names)
		obs[f'{coord}_0']=obs.index.to_series().map(df_coord[f'{coord}_0'].to_dict())
		obs[f'{coord}_1']=obs.index.to_series().map(df_coord[f'{coord}_1'].to_dict())
	print(obs.shape)
	obs.to_csv(os.path.join(workdir,"obs.tsv"),sep='\t')

	for level in ['Group','Subclass']:
		D=obs.groupby(level).apply(lambda x:str(x['anatomical_region'].value_counts().to_dict())).to_dict()
		if level=='Group':
			cols=['Neighborhood','Class','Subclass','Group']
		else:
			cols = ['Neighborhood', 'Class', 'Subclass']
		obs1=obs.reset_index().loc[:,cols].drop_duplicates()
		obs1['RegionDistrubution']=obs1[level].map(D)
		obs1.sort_values(cols,inplace=True)
		obs1.to_excel(os.path.join(workdir, f"{level}.taxonomy.xlsx"),index=False)

def mpl_style():
	import matplotlib as mpl
	mpl.style.use('default')
	mpl.rcParams['pdf.fonttype'] = 42
	mpl.rcParams['ps.fonttype'] = 42
	mpl.rcParams['font.family'] = 'sans-serif'
	mpl.rcParams['font.sans-serif'] = 'Arial'
	mpl.rcParams['figure.dpi'] = 80
	mpl.rcParams['savefig.dpi'] = 300

def parse_json(url):
	import json
	import requests
	data=json.loads(requests.get(url).content.decode())
	R=[]
	for record in data['msg']:
		R.append([record['acronym'],record['color_hex_triplet'],record['id'],record['name'],record['parent_structure_id'],record['safe_name'],record['structure_id_path']])
	df=pd.DataFrame(R,columns=['acronym','Hex','id','name','parent_structure_id','safe_name','structure_id_path'])
	return df

def get_brain_region_structure():
	"""
	https://atlas.brain-map.org/
	BICAN: https://atlas.brain-map.org/atlasviewer/ontologies/11.json
	HBA: https://atlas.brain-map.org/atlasviewer/ontologies/7.json
	Returns
	-------

	"""
	# pip install XlsxWriter
	writer = pd.ExcelWriter("AllenBrainRegionStructure.xlsx")
	for url,name in zip(['https://atlas.brain-map.org/atlasviewer/ontologies/11.json','https://atlas.brain-map.org/atlasviewer/ontologies/7.json'],['BICAN_Brodmann','HBA_guide']):
		print(name)
		df=parse_json(url)
		df['Hex']='#'+df.Hex.map(str)
		df['id']=df['id'].fillna('None').map(str)
		id2acronym=df.set_index('id').acronym.to_dict()
		df.insert(0,'Parent',df.parent_structure_id.fillna(0).map(int).map(str).map(id2acronym))
		df.parent_structure_id=df.parent_structure_id.map(str)
		df.structure_id_path=df.structure_id_path.apply(lambda x:x[1:-1])
		df['structure_path']=df.structure_id_path.apply(lambda x:'//'.join([id2acronym[p] for p in x.split('/') if p!='']))
		# data.style.background_gradient(cmap='gray_r')
		# data.style.applymap(lambda x:'color:'+x if x.startswith('#') else 'color: white')
		df=df.loc[:,['Parent','acronym','Hex','name','safe_name','id','parent_structure_id','structure_path','structure_id_path']]
		df.to_excel(writer,sheet_name=name,index=False)
		workbook = writer.book
		worksheet = writer.sheets[name]
		colors=df.Hex.tolist()
		col_idx=df.columns.tolist().index('Hex')
		for i in range(df.shape[0]):
			color = colors[i]
			f = workbook.add_format({'bold': True, 'font_color': 'black', 'bg_color': color})
			worksheet.write(i+1, col_idx, color,f) #worksheet.write(row, col, *args), row and col are 0-based
		# width=20
		# cell_fmt = workbook.add_format(
		#                 {'bold': False,'font_color': 'black',
		#                 # 'bg_color':'green',
		#                 'align': 'center', 'valign': 'vcenter'})
		# # styled = data.style.applymap(lambda val: 'color: %s' % 'red' if val < 0 else 'black').highlight_max()
		# worksheet.set_column(0,col_idx,width,cell_fmt) #first_col,last_col
		# worksheet.conditional_format(f'A:{last_col}', {'type': 'no_blanks', 'format': cell_fmt})
	writer.close()
	df_bican=pd.read_excel("Jon.xlsx")
	df_bican.Acronym=df_bican.Acronym.apply(lambda x:x.split('(')[0].strip())
	regions=df_bican.Acronym.tolist()
	df=pd.read_excel("AllenBrainRegionStructure.xlsx",
					 sheet_name="BICAN_Brodmann")
	# df['Keep']=df.structure_path.apply(lambda x:[r for r in x.split('//') if r in regions])
	# df=df.loc[df.Keep.apply(len) > 0]
	# df.drop('Keep',axis=1,inplace=True)
	df=df.loc[df.acronym.isin(regions)]
	df.parent_structure_id=df.parent_structure_id.map(int)

	writer = pd.ExcelWriter("BICAN_regions.xlsx")
	df.to_excel(writer, sheet_name='BICAN', index=False)
	workbook = writer.book
	worksheet = writer.sheets['BICAN']
	colors=df.Hex.tolist()
	col_idx=df.columns.tolist().index('Hex')
	for i in range(df.shape[0]):
		color = colors[i]
		f = workbook.add_format({'bold': True, 'font_color': 'black', 'bg_color': color})
		worksheet.write(i+1, col_idx, color,f)
	writer.close()

def read_google_sheet(url=None,**kwargs):
	assert not url is None
	# url="https://docs.google.com/spreadsheets/d/12H3p2F_qrcQ3ymF614VRzU_6vcVTsRca0uOIBQJXVaU/edit?gid=1969763406#gid=1969763406"
	Id=url.split('/d/')[1].split('/')[0]
	gid=url.split('?gid=')[1].split('#gid=')[0]
	df = pd.read_csv(f"https://docs.google.com/spreadsheets/d/{Id}/export?format=tsv&id={Id}&gid={gid}",
					 sep='\t',**kwargs)
	return df

def get_regions_mapping(project='BG',url=None):
	project_url=dict(BG="https://docs.google.com/spreadsheets/d/12H3p2F_qrcQ3ymF614VRzU_6vcVTsRca0uOIBQJXVaU/export?format=tsv&id=12H3p2F_qrcQ3ymF614VRzU_6vcVTsRca0uOIBQJXVaU&gid=0",
				    HIP="https://docs.google.com/spreadsheets/d/12H3p2F_qrcQ3ymF614VRzU_6vcVTsRca0uOIBQJXVaU/export?format=tsv&id=12H3p2F_qrcQ3ymF614VRzU_6vcVTsRca0uOIBQJXVaU&gid=2114739630",
					THM="https://docs.google.com/spreadsheets/d/12H3p2F_qrcQ3ymF614VRzU_6vcVTsRca0uOIBQJXVaU/export?format=tsv&id=12H3p2F_qrcQ3ymF614VRzU_6vcVTsRca0uOIBQJXVaU&gid=1291838884")
	if url is None:
		url=project_url[project]
	print(url)
	df = pd.read_csv(url, sep='\t')
	df.BICAN=df.BICAN.apply(lambda x:x.split(' ')[0] if not pd.isna(x) else np.nan)
	return df

def readTree(infile, outdir=None, name="BICAN", plot=True):
	import treelib
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	print("**-Reading tree")
	df = pd.read_excel(infile, sheet_name=name)
	tree = treelib.Tree()
	for i, row in df.iterrows():
		D = row.to_dict()
		branch = D['structure_path'].split('//')
		for j in range(len(branch)):
			Id = '//'.join(branch[:j + 1])
			if tree.contains(Id):
				continue
			data = dict(name=D['name'], fillcolor=D['Hex'], safe_name=D['safe_name'])
			if j == 0:
				node = tree.create_node(tag=branch[j], identifier=Id, data=data)
			else:
				parent_id='//'.join(branch[:j])
				node = tree.create_node(tag=branch[j], identifier=Id, parent=parent_id, data=data)

	tree.save2file(os.path.join(outdir, f'{name}.txt'))
	with open(os.path.join(outdir, f'{name}.json'), 'w', encoding='utf-8') as f:
		f.write(tree.to_json(with_data=True))
	with open(os.path.join(outdir, 'tree.pickle'), 'wb') as f:
		pickle.dump(tree, f, True)
	if plot:
		graphvizPlotTree(tree=tree, outdir=os.path.join(outdir, name))
	return tree

def tree2networkx(tree):
	import networkx as nx
	print("**-Converting tree to networkx")
	def get_edges(tree, root):
		if tree.depth(root) == tree.depth():
			return
		for child in tree.children(root.identifier):
			edges.append([root.identifier, child.identifier])  # identifier is better than tag.
			get_edges(tree, child)

	root = tree.get_node(tree.root)
	edges = []
	get_edges(tree, root)

	G = nx.DiGraph()
	G.add_edges_from(edges)
	for node in G.nodes():
		data=tree.get_node(node).data
		for k in data:
			G.nodes[node][k]=data[k]

	return G

def importPygraphviz():
	try:
		import pygraphviz
		from networkx.drawing.nx_agraph import graphviz_layout
	except ImportError:
		try:
			import pydotplus
			from networkx.drawing.nx_pydot import graphviz_layout
		except ImportError:
			raise ImportError("This function needs Graphviz and either \
							  PyGraphviz or PyDotPlus, Please install it \
							  from http://www.graphviz.org/ (window) or yum install graphviz / conda install pygraphviz(centos),\
							  then pip install pydotplus;pip install pydot.")
	return graphviz_layout

def plotNetworkx(G, outfile='network.pdf',prog='dot',figsize=(5, 6)):
	print("**-Plotting networkx")
	import networkx as nx
	from pym3c.utils import mpl_style
	mpl_style()
	# plt.switch_backend('agg')
	bin_path = os.path.dirname(sys.executable)
	os.environ['PATH'] = bin_path+':'+os.environ['PATH']
	graphviz_layout = importPygraphviz()

	fig = plt.figure(figsize=figsize)
	ax = plt.gca()
	pos = graphviz_layout(G, prog=prog)  # graphviz_layout(G, prog="twopi", root=0)
	#    pos=nx.spring_layout(G)
	degrees = [G.degree(n) for n in G.nodes()]
	d_max = max(degrees)
	d_min = min(degrees)
	d_gap = d_max - d_min

	node_sizes = [(3 + (d - d_min) / d_gap) * 10 for d in degrees]
	nx.draw_networkx_nodes(
		G, pos, nodelist=G.nodes(), node_size=node_sizes, \
		linewidths=0.1, vmin=0, vmax=1, alpha=0.5,node_shape='o',
		node_color=[node[1]['fillcolor'] for node in G.nodes(data=True)], ax=ax
	)

	nx.draw_networkx_labels(G, pos=pos,font_size=3,
							labels={node: node.split('//')[-1] for node in G.nodes()},
							horizontalalignment='center',ax=ax)
	nx.draw_networkx_edges(G, pos, edgelist=G.edges(), width=0.1, edge_color="black",
						   alpha=0.6, ax=ax, style='-', arrowsize=1, arrowstyle='-')

	#    plt.tight_layout()
	plt.savefig(outfile)
	plt.show()

def get_slab_mapping(
	url="https://docs.google.com/spreadsheets/d/1hOqa-bD2IlgdKuCxKq0UTdDSyCLPIfMfQ6Ab59mAXy0/edit?gid=1782718682"
):
	df=read_google_sheet(url,skiprows=1)
	df=df.iloc[:,3:].set_index('Distance').stack().reset_index()
	df.columns=['Distance','Donor','Slab']
	df=df.loc[df.Slab.notna()]
	df.Slab=df.Slab.apply(lambda x: x.split(' ')[0])
	D=df.groupby('Donor').apply(lambda x: x.set_index('Slab').Distance.to_dict()).to_dict()
	return D

def slab2distance(x,donor):
	"""
	Convert slab to distance
	:param x: slab name
	:return: distance
	"""
	import numpy as np
	if pd.isna(x) or x == '':
		return np.nan
	D= get_slab_mapping()
	if donor not in D:
		return np.nan
	distances=[]
	i=0
	prefix='CX'
	while i <= len(x)-1:
		e= x[i:i+2]
		if e in ['CX','BS','CB']:
			prefix=e
			i+=2
			continue
		slab_id= prefix + e
		if slab_id in D[donor]:
			distances.append(D[donor][slab_id])
		i+=2
	return np.mean(distances)

def despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False):
	"""
	Remove the top and right spines from plot(s).

	Parameters
	----------
	fig : matplotlib figure, optional
		Figure to despine all axes of, defaults to the current figure.
	ax : matplotlib axes, optional
		Specific axes object to despine. Ignored if fig is provided.
	top, right, left, bottom : boolean, optional
		If True, remove that spine.

	Returns
	-------
	None

	"""
	if fig is None and ax is None:
		axes = plt.gcf().axes
	elif fig is not None:
		axes = fig.axes
	elif ax is not None:
		axes = [ax]

	for ax_i in axes:
		for side in ["top", "right", "left", "bottom"]:
			is_visible = not locals()[side]
			ax_i.spines[side].set_visible(is_visible)
		if left and not right:  # remove left, keep right
			maj_on = any(t.tick1line.get_visible() for t in ax_i.yaxis.majorTicks)
			min_on = any(t.tick1line.get_visible() for t in ax_i.yaxis.minorTicks)
			ax_i.yaxis.set_ticks_position("right")
			for t in ax_i.yaxis.majorTicks:
				t.tick2line.set_visible(maj_on)
			for t in ax_i.yaxis.minorTicks:
				t.tick2line.set_visible(min_on)

		if bottom and not top:
			maj_on = any(t.tick1line.get_visible() for t in ax_i.xaxis.majorTicks)
			min_on = any(t.tick1line.get_visible() for t in ax_i.xaxis.minorTicks)
			ax_i.xaxis.set_ticks_position("top")
			for t in ax_i.xaxis.majorTicks:
				t.tick2line.set_visible(maj_on)
			for t in ax_i.xaxis.minorTicks:
				t.tick2line.set_visible(min_on)

def _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=5):
    # This function assume coord is [0, 1].
    # clean ax axises
    ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    despine(ax=ax, left=True, bottom=True)

    _arrow_kws = {"width": 0.003, "linewidth": 0, "color": "black"}
    if arrow_kws is not None:
        _arrow_kws.update(arrow_kws)

    ax.arrow(0.06, 0.06, 0, 0.06, **_arrow_kws, transform=ax.transAxes)
    ax.arrow(0.06, 0.06, 0.06, 0, **_arrow_kws, transform=ax.transAxes)
    ax.text(
        0.06,
        0.03,
        x.upper().replace("_", " "),
        fontdict={"fontsize": fontsize, "horizontalalignment": "left", "verticalalignment": "center"},
        transform=ax.transAxes,
    )
    ax.text(
        0.03,
        0.06,
        y.upper().replace("_", " "),
        fontdict={
            "fontsize": fontsize,
            "rotation": 90,
            "rotation_mode": "anchor",
            "horizontalalignment": "left",
            "verticalalignment": "center",
        },
        transform=ax.transAxes,
    )
    return

def zoom_min_max(vmin, vmax, scale):
    """Zoom min and max value."""
    width = vmax - vmin
    width_zoomed = width * scale
    delta_value = (width_zoomed - width) / 2
    return vmin - delta_value, vmax + delta_value

def zoom_ax(ax, zoom_scale, on="both"):
    """Zoom ax on both x and y-axis."""
    on = on.lower()
    xlim = ax.get_xlim()
    xlim_zoomed = zoom_min_max(vmin=xlim[0], vmax=xlim[1],scale=zoom_scale)

    ylim = ax.get_ylim()
    ylim_zoomed = zoom_min_max(vmin=ylim[0], vmax=ylim[1],scale=zoom_scale)

    if (on == "both") or ("x" in on):
        ax.set_xlim(xlim_zoomed)
    if (on == "both") or ("y" in on):
        ax.set_ylim(ylim_zoomed)

def _extract_coords(data, coord_base, x, y):
	import xarray as xr
	import anndata
	if (x is not None) and (y is not None):
		pass
	else:
		x = f"{coord_base}_0"
		y = f"{coord_base}_1"

	if isinstance(data, anndata.AnnData):
		adata = data
		_data = pd.DataFrame(
			{
				"x": adata.obsm[f"X_{coord_base}"][:, 0],
				"y": adata.obsm[f"X_{coord_base}"][:, 1],
			},
			index=adata.obs_names,
		)
	elif isinstance(data, xr.Dataset):
		ds = data
		if coord_base not in ds.dims:
			raise KeyError(f"xr.Dataset do not contain {coord_base} dim")
		data_var = {i for i in ds.data_vars.keys() if i.startswith(coord_base)}.pop()
		_data = pd.DataFrame(
			{
				"x": ds[data_var].sel({coord_base: f"{coord_base}_0"}).to_pandas(),
				"y": ds[data_var].sel({coord_base: f"{coord_base}_1"}).to_pandas(),
			}
		)
	else:
		if (x not in data.columns) or (y not in data.columns):
			raise KeyError(f"{x} or {y} not found in columns.")
		_data = pd.DataFrame({"x": data[x], "y": data[y]})
	return _data, x, y

def _density_based_sample(data: pd.DataFrame, coords: list, portion=None, size=None, seed=None):
	"""Down sample data based on density, to prevent overplot in dense region and decrease plotting time."""
	from sklearn.neighbors import LocalOutlierFactor
	clf = LocalOutlierFactor(
		n_neighbors=20,
		algorithm="auto",
		leaf_size=30,
		metric="minkowski",
		p=2,
		metric_params=None,
		contamination=0.1,
	)

	# coords should already exist in data, get them by column names list
	data_coords = data[coords]
	clf.fit(data_coords)
	# original score is negative, the larger the denser
	density_score = clf.negative_outlier_factor_
	delta = density_score.max() - density_score.min()
	# density score to probability: the denser the less probability to be picked up
	probability_score = 1 - (density_score - density_score.min()) / delta
	probability_score = np.sqrt(probability_score)
	probability_score = probability_score / probability_score.sum()

	if size is not None:
		pass
	elif portion is not None:
		size = int(data_coords.index.size * portion)
	else:
		raise ValueError("Either portion or size should be provided.")
	if seed is not None:
		np.random.seed(seed)
	selected_cell_index = np.random.choice(
		data_coords.index, size=size, replace=False, p=probability_score
	)  # choice data based on density weights

	# return the down sampled data
	return data.reindex(selected_cell_index)

def _auto_size(ax, n_dots):
    """Auto determine dot size based on ax size and n dots"""
    bbox = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted())
    scale = bbox.width * bbox.height / 14.6  # 14.6 is a 5*5 fig I used to estimate
    n = n_dots / scale  # larger figure means data look sparser
    if n < 500:
        s = 14 - n / 100
    elif n < 1500:
        s = 7
    elif n < 3000:
        s = 5
    elif n < 8000:
        s = 3
    elif n < 15000:
        s = 2
    elif n < 30000:
        s = 1.5
    elif n < 50000:
        s = 1
    elif n < 80000:
        s = 0.8
    elif n < 150000:
        s = 0.6
    elif n < 300000:
        s = 0.5
    elif n < 500000:
        s = 0.4
    elif n < 800000:
        s = 0.3
    elif n < 1000000:
        s = 0.2
    elif n < 2000000:
        s = 0.1
    elif n < 3000000:
        s = 0.07
    elif n < 4000000:
        s = 0.05
    elif n < 5000000:
        s = 0.03
    else:
        s = 0.02
    return s

def _take_data_series(data, k):
	import xarray as xr
	import anndata
	if isinstance(data, (xr.Dataset, xr.DataArray)):
		_value = data[k].to_pandas()
	elif isinstance(data, anndata.AnnData):
		_value = data.obs[k].copy()
	else:
		_value = data[k].copy()
	return _value

def level_one_palette(name_list, order=None, palette="auto"):
	name_set = set(name_list.dropna())
	if palette == "auto":
		if len(name_set) < 10:
			palette = "tab10"
		elif len(name_set) < 20:
			palette = "tab20"
		else:
			palette = "rainbow"

	if order is None:
		try:
			order = sorted(name_set)
		except TypeError:
			# name set contains multiple dtype (e.g., str and np.NaN)
			# do not sort
			order = list(name_set)
	else:
		if (set(order) != name_set) or (len(order) != len(name_set)):
			raise ValueError("Order is not equal to set(name_list).")

	n = len(order)
	colors = sns.color_palette(palette, n)
	color_palette = {}
	for name, color in zip(order, colors):
		color_palette[name] = color
	return color_palette

def _calculate_luminance(color):
	"""
	Calculate the relative luminance of a color according to W3C standards

	Parameters
	----------
	color : matplotlib color or sequence of matplotlib colors
		Hex code, rgb-tuple, or html color name.
	Returns
	-------
	luminance : float(s) between 0 and 1

	"""
	rgb = matplotlib.colors.colorConverter.to_rgba_array(color)[:, :3]
	rgb = np.where(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
	lum = rgb.dot([0.2126, 0.7152, 0.0722])
	try:
		return lum.item()
	except ValueError:
		return lum

def _text_anno_scatter(
    data: pd.DataFrame,
    ax,
    x: str,
    y: str,
    dodge_text=False,
    anno_col="text_anno",
    text_kws=None,
    text_transform=None,
    dodge_kws=None,
    luminance=0.48
):
	"""Add text annotation to a scatter plot."""
	import copy
	# prepare kws
	text_kws={} if text_kws is None else text_kws
	text_kws.setdefault("fontsize",5)
	text_kws.setdefault("fontweight","black")
	text_kws.setdefault("ha","center") #horizontalalignment
	text_kws.setdefault("va","center") #verticalalignment
	text_kws.setdefault("color","black") #c
	bbox=dict(boxstyle='round',edgecolor=(0.5, 0.5, 0.5, 0.2),fill=False,
								facecolor=(0.8, 0.8, 0.8, 0.2),alpha=1,linewidth=0.5)
	text_kws.setdefault("bbox",bbox)
	for key in bbox:
		if key not in text_kws['bbox']:
			text_kws['bbox'][key]=bbox[key]
	# plot each text
	text_list = []
	for text, sub_df in data.groupby(anno_col):
		if text_transform is None:
			text = str(text)
		else:
			text = text_transform(text)
		if text.lower() in ["", "nan"]:
			continue
		_x, _y = sub_df[[x, y]].median()
		
		use_text_kws=copy.deepcopy(text_kws) #text_kws.copy()
		if isinstance(text_kws['bbox']['facecolor'],dict):
			use_text_kws['bbox']['facecolor']=text_kws['bbox']['facecolor'].get(text,'gray')
		if isinstance(text_kws['color'],dict):
			use_color=text_kws['color'].get(text,'black')
			use_text_kws['color']=use_color
		if not luminance is None and not use_text_kws['bbox']['facecolor'] is None:
			lum = _calculate_luminance(use_text_kws['bbox']['facecolor'])
			if lum > luminance:
				use_text_kws['color']='black'
				use_text_kws['bbox']['edgecolor']='black'
		
		text = ax.text(
			_x,
			_y,
			text,
			**use_text_kws
		)
		text_list.append(text)

	if dodge_text:
		try:
			from adjustText import adjust_text

			_dodge_parms = {
				"force_points": (0.02, 0.05),
				"arrowprops": {
					"arrowstyle": "->",
					"fc": "black",
					"ec": "none",
					"connectionstyle": "angle,angleA=-90,angleB=180,rad=5",
				},
				"autoalign": "xy",
			}
			if dodge_kws is not None:
				_dodge_parms.update(dodge_kws)
			adjust_text(text_list, x=data["x"], y=data["y"], **_dodge_parms)
		except ModuleNotFoundError:
			print("Install adjustText package to dodge text, see its github page for help")
	return

def tight_hue_range(hue_data, portion):
	"""Automatic select a SMALLEST data range that covers [portion] of the data."""
	hue_data = hue_data[np.isfinite(hue_data)]
	hue_quantiles = hue_data.quantile(q=np.arange(0, 1, 0.01))
	min_window_right = (
		hue_quantiles.rolling(window=int(portion * 100)).apply(lambda i: i.max() - i.min(), raw=True).idxmin()
	)
	min_window_left = max(0, min_window_right - portion)
	vmin, vmax = tuple(hue_data.quantile(q=[min_window_left, min_window_right]))
	if np.isfinite(vmin):
		vmin = max(hue_data.min(), vmin)
	else:
		vmin = hue_data.min()
	if np.isfinite(vmax):
		vmax = min(hue_data.max(), vmax)
	else:
		vmax = hue_data.max()
	return vmin, vmax

def density_contour(
    ax,
    data,
    x,
    y,
    groupby=None,
    c="lightgray",
    single_contour_pad=1,
    linewidth=1,
    palette=None,
):
	from sklearn.neighbors import LocalOutlierFactor
	_data = data.copy()

	if groupby is not None:
		if isinstance(groupby, str):
			_data["groupby"] = data[groupby]
		else:
			_data["groupby"] = groupby
	else:
		_data["groupby"] = "one group"

	_contour_kws = {"linewidths": linewidth, "levels": (-single_contour_pad,), "linestyles": "dashed"}
	_lof_kws = {"n_neighbors": 25, "novelty": True, "contamination": "auto"}

	xmin, ymin = _data[[x, y]].min()
	xmax, ymax = _data[[x, y]].max()
	xmin, xmax = zoom_min_max(xmin, xmax, 1.2)
	ymin, ymax = zoom_min_max(ymin, ymax, 1.2)

	for group, sub_data in _data[[x, y, "groupby"]].groupby("groupby"):
		xx, yy = np.meshgrid(np.linspace(xmin, xmax, 500), np.linspace(ymin, ymax, 500))
		clf = LocalOutlierFactor(**_lof_kws)
		clf.fit(sub_data.iloc[:, :2].values)
		z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
		z = z.reshape(xx.shape)
		if palette is None:
			_color = c
		else:
			_color = palette[group] if group in palette else c
		# plot contour line(s)
		ax.contour(xx, yy, z, colors=_color, **_contour_kws)
	return

def plot_color_dict_legend(
	D, ax=None, title=None, color_text=True, 
	kws=None,luminance=0.5
):
	"""
	plot legned for color dict

	Parameters
	----------
	D: a dict, key is categorical variable, values are colors.
	ax: axes to plot the legend.
	title: title of legend.
	color_text: whether to change the color of text based on the color in D.
	label_side: right of left.
	kws: kws passed to plt.legend.

	Returns
	-------
	ax.legend

	"""
	import matplotlib.patches as mpatches
	if ax is None:
		ax = plt.gca()
	lgd_kws = kws.copy() if not kws is None else {}  # bbox_to_anchor=(x,-0.05)
	lgd_kws.setdefault("frameon", True)
	lgd_kws.setdefault("ncol", 1)
	lgd_kws["loc"] = "upper left"
	lgd_kws.setdefault("borderpad", 0.1 * mm2inch * 72)  # 0.1mm
	lgd_kws.setdefault("markerscale", 1)
	lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
	lgd_kws.setdefault("handlelength", 1)  # font size, units is points
	lgd_kws.setdefault(
		"borderaxespad", 0.1
	)  # The pad between the axes and legend border, in font-size units.
	lgd_kws.setdefault(
		"handletextpad", 0.4
	)  # The pad between the legend handle and text, in font-size units.
	lgd_kws.setdefault(
		"labelspacing", 0.15
	)  # gap height between two Patches,  0.05*mm2inch*72
	lgd_kws.setdefault("columnspacing", 0.5)
	# lgd_kws["bbox_transform"] = ax.figure.transFigure
	lgd_kws.setdefault("bbox_to_anchor", (1, 1))
	lgd_kws.setdefault("title", title)
	lgd_kws.setdefault("markerfirst", True)
	l = [
		mpatches.Patch(color=c, label=l) for l, c in D.items()
	]  # kws:?mpatches.Patch; rasterized=True
	ms = lgd_kws.pop("markersize", 10)
	L = ax.legend(handles=l, **lgd_kws)
	L._legend_box.align = 'center'
	L.get_title().set_ha('center')
	if color_text:
		for text in L.get_texts():
			try:
				lum = _calculate_luminance(D[text.get_text()])
				if luminance is None:
					text_color = "black"
				else:
					text_color = "black" if lum > luminance else D[text.get_text()]
				text.set_color(text_color)
			except:
				pass
	# ax.add_artist(L)
	ax.grid(False)
	return L

def plot_marker_legend(
	color_dict=None, ax=None, title=None, color_text=True, 
	marker='o',kws=None,luminance=0.5
):
	"""
	plot legned for different marker

	Parameters
	----------
	D: a dict, key is categorical variable, values are marker.
	ax: axes to plot the legend.
	title: title of legend.
	color_text: whether to change the color of text based on the color in D.
	label_side: right of left.
	kws: kws passed to plt.legend.

	Returns
	-------
	ax.legend

	"""
	if ax is None:
		ax = plt.gca()

	lgd_kws = kws.copy() if not kws is None else {}  # bbox_to_anchor=(x,-0.05)
	lgd_kws.setdefault("frameon", True)
	lgd_kws.setdefault("ncol", 1)
	lgd_kws["loc"] = "upper left"
	# lgd_kws["bbox_transform"] = ax.figure.transFigure
	lgd_kws.setdefault("borderpad", 0.2 * mm2inch * 72)  # 0.1mm
	# lgd_kws.setdefault("markerscale", 1)
	lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
	lgd_kws.setdefault("handlelength", 1)  # font size, units is points
	lgd_kws.setdefault(
		"borderaxespad", 0.1
	)  # The pad between the axes and legend border, in font-size units.
	lgd_kws.setdefault(
		"handletextpad", 0.4 #0.2 * mm2inch * 72
	)  # The pad between the legend handle and text, in font-size units.
	lgd_kws.setdefault(
		"labelspacing", 0.15
	)  # gap height between two Patches,  0.05*mm2inch*72
	lgd_kws.setdefault("columnspacing", 0.5)
	lgd_kws.setdefault("bbox_to_anchor", (1, 1))
	lgd_kws.setdefault("title", title)
	lgd_kws.setdefault("markerfirst", True)

	ms = lgd_kws.pop("markersize", 10)
	L = [
		Line2D(
			[0],
			[0],
			color=color,
			marker=marker,
			linestyle="None",
			markersize=ms,
			label=l,
		)
		for l,color in color_dict.items()
	]
	L = ax.legend(handles=L, **lgd_kws)
	ax.figure.canvas.draw()
	L._legend_box.align = 'center'
	L.get_title().set_ha('center')
	if color_text:
		for text in L.get_texts():
			try:
				lum = _calculate_luminance(color_dict[text.get_text()])
				if luminance is None:
					text_color = "black"
				else:
					text_color = "black" if lum > luminance else color_dict[text.get_text()]
				text.set_color(text_color)
			except:
				pass
	# ax.add_artist(lgd)
	ax.grid(False)
	return L

# Custom handler for legend: circle text as marker + label
class TextWithCircleHandler(HandlerBase):
	def __init__(self, marker_text='', label_text='', 
			  text_kws={}, **kwargs):
		HandlerBase.__init__(self, **kwargs)
		self.marker_text = marker_text
		self.text_kws=text_kws

	def create_artists(self, legend, orig_handle,
						xdescent, ydescent, width, height, fontsize, trans):
		# Marker (number with circle)
		self.text_kws.setdefault("fontsize",fontsize)
		# print(self.text_kws)
		shift=2 * self.text_kws['fontsize'] * 0.65 / 72 / mm2inch
		circ_text = Text(
			xdescent + legend.borderaxespad + shift,  height / 2,
			self.marker_text, 
			**self.text_kws
		)
		return [circ_text]
	
def plot_text_legend(color_dict, code2label, ax=None, title=None, color_text=True, 
					 boxstyle='Circle',marker_pad=0.1,legend_kws=None,marker_fontsize=4,
					 text_kws=None,alpha=0.7,luminance=0.5):
	import copy
	# print(color_dict)
	lgd_kws = legend_kws.copy() if not legend_kws is None else {}  # bbox_to_anchor=(x,-0.05)
	lgd_kws.setdefault("frameon", True)
	lgd_kws.setdefault("ncol", 1)
	lgd_kws["loc"] = "upper left"
	# lgd_kws["bbox_transform"] = ax.figure.transFigure
	lgd_kws.setdefault("borderpad", 0.1 * mm2inch * 72)  # 0.1mm
	# lgd_kws.setdefault("markerscale", 1)
	lgd_kws.setdefault("handleheight", 0.5)  # font size, units is points
	lgd_kws.setdefault("handlelength", 1)  # font size, units is points
	lgd_kws.setdefault(
		"borderaxespad", 0.5
	)  # The pad between the axes and legend border, in font-size units.
	lgd_kws.setdefault(
		"handletextpad", 0.4
	)  # The pad between the legend handle and text, in font-size units.
	lgd_kws.setdefault(
		"labelspacing", 0.2
	)  # gap height between two row of legend,  0.05*mm2inch*72
	lgd_kws.setdefault("columnspacing", 0.5)
	lgd_kws.setdefault("bbox_to_anchor", (1, 1))
	lgd_kws.setdefault("title", title)
	lgd_kws.setdefault("markerfirst", True)
	ms = lgd_kws.pop("markersize", 10)

	# text_kws
	if text_kws is None:
		text_kws={}
	default_marker_text_kws=dict(
		bbox=dict(boxstyle=f"{boxstyle},pad={marker_pad}", #Square, Circle, Round
			edgecolor='black',linewidth=0.4,
			fill=True,facecolor='white',alpha=alpha),
		horizontalalignment='center', verticalalignment='center',
		fontsize=marker_fontsize,color='black')
	for k in default_marker_text_kws:
		if k == 'bbox':
			if k not in text_kws:
				text_kws['bbox']=default_marker_text_kws[k]
			else:
				for k1 in default_marker_text_kws['bbox'].keys():
					if k1 not in text_kws['bbox']:
						text_kws['bbox'].setdefault(k1,default_marker_text_kws['bbox'][k1])
		if k not in text_kws:
			text_kws.setdefault(k,default_marker_text_kws[k])
	
	# Create handles and handlers
	handles = []
	handler_map = {}
	for code in sorted([int(i) for i in code2label.keys()]):
		label=code2label[code]
		code_text=str(code)
		handle = Line2D([], [], linestyle=None,
				  label=label)
		handles.append(handle)
		color=color_dict.get(label, 'black')
		text_kws1= copy.deepcopy(text_kws)
		text_kws1['bbox']['facecolor']=color
		lum = _calculate_luminance(color)
		if lum <= 0.1: # for black-like color, use white marker text
			text_kws1['color']='white'
		# print(bbox)
		handler_map[handle] = TextWithCircleHandler(
			marker_text=code_text, label_text=label,
			text_kws=text_kws1,
			)

	# Draw custom legend
	L=ax.legend(handles=handles, handler_map=handler_map, 
			 **lgd_kws)
	L._legend_box.align = 'center'
	L.get_title().set_ha('center')
	ax.figure.canvas.draw()
	if color_text:
		for text in L.get_texts():
			try:
				lum = _calculate_luminance(color_dict[text.get_text()])
				if luminance is None:
					text_color = "black"
				else:
					text_color = "black" if lum > luminance else color_dict[text.get_text()]
				text.set_color(text_color)
			except:
				pass
	# ax.add_artist(lgd)
	# ax.grid(False)

def plot_cmap_legend(
	cax=None, ax=None, cmap="turbo", label=None, kws=None,
	labelsize=6, linewidth=0.5,ticklabel_size=4,
):
	"""
	Plot legend for cmap.

	Parameters
	----------
	cax : Axes into which the colorbar will be drawn.
	ax :  axes to anchor.
	cmap : turbo, hsv, Set1, Dark2, Paired, Accent,tab20,exp1,exp2,meth1,meth2
	label : title for legend.
	kws : dict
		kws passed to plt.colorbar (matplotlib.figure.Figure.colorbar).

	Returns
	-------
	cbar: axes of legend

	"""
	label = "" if label is None else label
	cbar_kws = {} if kws is None else kws.copy()
	cbar_kws.setdefault("label", label)
	# cbar_kws.setdefault("aspect",3)
	cbar_kws.setdefault("orientation", "vertical")
	# cbar_kws.setdefault("use_gridspec", True)
	# cbar_kws.setdefault("location", "bottom")
	cbar_kws.setdefault("fraction", 1)
	cbar_kws.setdefault("shrink", 1)
	cbar_kws.setdefault("pad", 0)
	cbar_kws.setdefault("extend", 'both')
	cbar_kws.setdefault("extendfrac", 0.1)
	# print(cbar_kws,kws)
	# print(type(cax))
	vmax = cbar_kws.pop("vmax", 1)
	vmin = cbar_kws.pop("vmin", 0)
	# print(vmin,vmax,'vmax,vmin')
	cax.set_ylim([vmin, vmax])
	# print(cax.get_ylim())
	vcenter= (vmax + vmin) / 2
	center=cbar_kws.pop("center",None)
	if center is None:
		center=vcenter
		m = plt.cm.ScalarMappable(
			norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap
		)
	else:
		m = plt.cm.ScalarMappable(
			norm=matplotlib.colors.TwoSlopeNorm(center,vmin=vmin, vmax=vmax), cmap=cmap
		)
	cbar_kws.setdefault("ticks", [vmin, center, vmax])
	cax.yaxis.set_label_position('right')
	cax.yaxis.set_ticks_position('right')
	cbar = ax.figure.colorbar(m, cax=cax, **cbar_kws)  # use_gridspec=True
	cbar.ax.tick_params(labelsize=ticklabel_size, size=ticklabel_size,width=linewidth) # size is for ticks, labelsize is for the number on ticks (ticklabels)
	cbar.ax.yaxis.label.set_fontsize(labelsize) # colorbar title fontsize
	cbar.ax.grid(False)
	return cbar

def normalize_mc_by_cell(
		use_adata,normalize_per_cell=True,
		clip_norm_value=10,verbose=1,hypo_score=False):
	from scipy.sparse import issparse
	normalized_flag = use_adata.uns.get('normalize_per_cell',False)
	if normalize_per_cell and not normalized_flag:  # divide frac by prior mean (determined by alpha and beta) for each cell
		# get normalized X
		cols = use_adata.obs.columns.tolist()
		if 'prior_mean' in cols:
			if verbose > 0:
				print("Normalizing cell level fraction by alpha and beta (prior_mean)")
			na_sum = use_adata.to_df().isna().sum().sum()
			if na_sum > 0:
				D = use_adata.obs.prior_mean.to_dict()
				if not hypo_score:
					use_adata.X = use_adata.to_df().apply(lambda x:x.fillna(D[x.name]) / D[x.name],axis=1).values
				else: # hypo score, the larger the value, the more hypomethylated
					use_adata.X = use_adata.to_df().apply(lambda x:D[x.name] / x.fillna(D[x.name]),axis=1).values
			else:
				if not hypo_score: # the smaller the value, the lower of the methylation fraction
					use_adata.X = use_adata.X / use_adata.obs.prior_mean.values[:, None]  # range = [0,1,10]
				else: # hypo-score
					use_adata.X = use_adata.obs.prior_mean.values[:, None] / use_adata.X
			if not clip_norm_value is None:
				if issparse(use_adata.X):
					X=use_adata.X.toarray()
				else:
					X=use_adata.X
				use_adata.X = np.clip(X, None, clip_norm_value)
			use_adata.uns['normalize_per_cell'] = True
		else:
			if verbose > 0:
				print("'prior_mean' not found in obs")
	elif normalize_per_cell and normalized_flag:
		logger.info("Input adata is already normalized, skip normalize_per_cell !")
	else:
		pass
	return use_adata
