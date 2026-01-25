# -*- coding: utf-8 -*-
"""
@author: DingWB
"""
import os
import pandas as pd
import numpy as np
import scanpy as sc
from .utils import serialize,normalize_mc_by_cell
from .tools import parse_gtf,cal_tpm
import anndata
from concurrent.futures import ProcessPoolExecutor, as_completed
from loguru import logger as logger

def merge_adata_regions(
		pseudobulk_adata_path,bin_size=5000,use_raw=True,
		res=100000,filter_chroms=True,boundary=None,
		exclude_chroms=['chrY','chrM','chrX']):
	"""
	Aggregrate 5kb RNA or ATAC adata into 100kb (25kb, or 10kb)

	Parameters
	----------
	pseudobulk_adata_path : path
	bin_size: int
	use_raw: 
	res : int
		100000
	filter_chroms : bool
		True

	Returns
		A dataframe
	-------

	"""
	import scanpy as sc
	def assign_boundary(indexes, names, list_x):
		"""
        For example, the 1st 25kb: 0,1,2,3,4; 2nd 25kb: 0,1,2,3,4; the boundary between these two bins should be: 3,4,0,1
        """
		name2index = {}
		boundary = []
		flag = True
		tmp_chrom = None
		tmp_name = ''
		for x, idx, name in zip(list_x, indexes, names):
			chrom = name.split('_')[0]
			if chrom != tmp_chrom:
				if len(boundary) > 0:
					name2index[tmp_name] = tuple(boundary)
				boundary = []
				flag = True
			if x == 0 or x == 1:
				boundary.append(idx)
				if x == 1:
					if flag:
						name2index[name] = tuple(boundary)
						flag = False
					if len(boundary) == 4:
						name2index[name] = tuple(boundary)
			if x == 2:
				boundary = []
				tmp_chrom = chrom
			if x == 3 or x == 4:
				boundary.append(idx)
			tmp_chrom = chrom
			tmp_name = name
		return name2index

	if isinstance(pseudobulk_adata_path,str):
		adata=anndata.read_h5ad(os.path.expanduser(pseudobulk_adata_path))
	else:
		adata=pseudobulk_adata_path
	if use_raw:
		if not adata.raw is None:
			adata=adata.raw.to_adata()
		else:
			logger.warning("adata.raw is None!!")
	data=adata.to_df().T
	fea="100kb" if res==100000 else '10kb' if res==10000 else '25kb'
	if boundary is None:
		boundary=True if fea == '25kb' else False
	groups=data.columns.tolist()
	data['chrom']=data.index.to_series().apply(lambda x:x.split(':')[0])
	data['start']=data.index.to_series().apply(lambda x:x.split(':')[1].split('-')[0]).map(int)
	data['BinID']=data['start'].apply(lambda x:np.floor(x / bin_size)).astype(int)

	if not boundary:
		# merge 5kb into 100kb or 10kb
		data[fea] = data.apply(lambda x:x['chrom']+'_'+str(x['start'] // res),axis=1)
		# index, for example, chr1_0, chr1_0, chr1_0...,chr1_1
		data = data.loc[:,groups+[fea]].groupby(fea).sum().T
	else: # for 25kb, get the domain doundary (10kb), not the 25kb windows.
		print("Generating results for boundaries of 25kb bins")
		ids = data['BinID'] % 5  # [0,1],2,[3,4,0,1],2,[3,4,0,1],...
		data[fea] = data.apply(lambda x:x['chrom']+'_'+str(x['BinID'] // 5),axis=1)
		name2index = assign_boundary(data.index.tolist(), data[fea].tolist(), ids.tolist())
		idx2name = {}
		for name in name2index:
			for idx in name2index[name]:
				idx2name[idx] = name
		data['NAME'] = data.index.to_series().map(idx2name)
		data = data.loc[data['NAME'].notna()]
		data = data.loc[data['NAME'].apply(lambda x: x.split('_')[0]) == data.chrom]
		data = data.loc[:,groups+['NAME']].groupby('NAME').sum().T
	if use_raw: # convert sum of raw counts into CPM
		adata = anndata.AnnData(X=data)
		sc.pp.normalize_total(adata, target_sum=1e6)
		sc.pp.log1p(adata) # log(CPM)
		data=adata.to_df()
	if filter_chroms:
		s_col=data.columns.to_frame()
		s_col['chrom']=s_col.index.to_series().apply(lambda x:x.split('_')[0])
		keep_cols=s_col.loc[s_col['chrom'].apply(lambda x:x not in exclude_chroms and len(x)<6)].index.tolist()
		# use_rows = list(set(mc_df.index.tolist()) & set(cellclass2majortype[group]))
		data = data.loc[:, keep_cols]
	# rows are 100kb bins, columns are cell types
	return data # to do next: subset rows (df_bin.index) and columns (cell types order)

def cal_stats(df_data,modality="RNA",
			  expression_cutoff=0):
	# Compute per-gene statistics (min, q25, q50, q75, max, sum) across cells.
	# Use NumPy's nanpercentile and nansum which are fast (uses quickselect under the hood).
	qs = np.nanpercentile(df_data.values, [0, 25, 50, 75, 100], axis=0)
	sums = np.nansum(df_data.values, axis=0) # for each column

	# fraction of cells expressing (or hypomethylated) the gene
	if modality!='RNA': # methylation, cutoff = 1
		# frac = df_data.apply(lambda x: x[x < 1].shape[0] / x.shape[0])
		# vectorized: count values < 1 per column divided by number of cells
		frac = (df_data < 1).sum(axis=0) / float(df_data.shape[0])
	else: # for RNA
		# frac = df_data.apply(lambda x: x[x > expression_cutoff].shape[0] / x.shape[0])
		# vectorized: count values > cutoff per column divided by number of cells
		frac = (df_data > expression_cutoff).sum(axis=0) / float(df_data.shape[0])
	return qs,sums,frac,df_data.columns.tolist()

def to_pseudobulk(
	adata_path,downsample=None,
	obs_path=None,groupby="Group",use_raw=True,expression_cutoff=0,
	modality="RNA",n_jobs=1,
	normalize_per_cell=True,clip_norm_value=10,
	normalization=None,target_sum=1e6,gtf=None,save=None
):
	if modality!='RNA': # methylation
		assert normalize_per_cell==True, "For methylation, normalize_per_cell should be True"
	raw_adata=anndata.read_h5ad(os.path.expanduser(adata_path),backed='r')
	if not obs_path is None:
		if isinstance(obs_path,str):
			obs=pd.read_csv(os.path.expanduser(obs_path),
				sep='\t',index_col=0)
		else:
			obs=obs_path.copy()
		overlapped_cells=list(set(raw_adata.obs_names.tolist()) & set(obs.index.tolist()))
		obs=obs.loc[overlapped_cells]
		raw_adata.obs[groupby]=raw_adata.obs.index.to_series().map(obs[groupby].to_dict())
	if not downsample is None:
		all_cells = raw_adata.obs.loc[raw_adata.obs[groupby].notna()].groupby(groupby).apply(
					lambda x: x.sample(downsample).index.tolist() if x.shape[0] > downsample else x.index.tolist()).sum()
	else:
		all_cells=raw_adata.obs.loc[raw_adata.obs[groupby].notna()].index.tolist()
	data={}
	if n_jobs==-1:
		n_jobs=os.cpu_count()
	with ProcessPoolExecutor(n_jobs) as executor:
		futures = {}
		for group in raw_adata.obs.loc[all_cells,groupby].unique():
			obs1=raw_adata.obs.loc[all_cells]
			use_cells=obs1.loc[obs1[groupby]==group].index.tolist()
			if len(use_cells)==0:
				continue
			adata=raw_adata[use_cells,:].to_memory()
			if modality!='RNA':
				adata = normalize_mc_by_cell(
					use_adata=adata, normalize_per_cell=normalize_per_cell,
					clip_norm_value=clip_norm_value,
					hypo_score=False,verbose=0)
			else:
				if use_raw and not adata.raw is None:
					adata.X=adata.raw.X.copy()
			df_data=adata.to_df() # rows are cells, columns are genes
			future = executor.submit(
				cal_stats,df_data=df_data,modality=modality,
				expression_cutoff=expression_cutoff
			)
			futures[future] = group
		logger.debug(f"Submitted {len(futures)} groups for pseudobulk calculation.")
		for future in as_completed(futures):
			group = futures[future]
			logger.debug(group)
			qs,sums,frac,header = future.result()
			for k,v in zip(['min', 'q25', 'q50', 'q75', 'max', 'sum'], qs.tolist() + [sums.tolist()]):
				if k not in data:
					data[k] = []
				data[k].append(pd.Series(v, name=group, index=header))
			frac.name = group
			if 'frac' not in data:
				data['frac'] = []
			data['frac'].append(pd.Series(frac, name=group,index=header))
	raw_adata.file.close()
	# for RNA, put sum of raw counts into adata.X
	X=pd.concat(data['sum'],axis=1).T # sum of raw counts or normalized methylation fraction
	vc=raw_adata.obs.loc[all_cells][groupby].value_counts().to_frame(name='cell_count')
	vc_dict=vc.to_dict()['cell_count']
	if modality!='RNA': # for methylation, put mean methylation level into adata.X
		X=X.apply(lambda x:x/vc_dict[x.name],axis=1)
	adata = anndata.AnnData(X=X,obs=vc.loc[X.index.tolist()])
	for k in data:
		if k=='sum':
			continue
		adata.layers[k]=pd.concat(objs=data[k],axis=1).T
	del data

	if not normalization is None:
		# Calculate CPM or TPM only if aggfunc is sum
		logger.info(f"Normalizing pseudobulk adata using {normalization} method.")
		if not gtf is None:
			df_gene = parse_gtf(gtf=gtf)
			# ['chrom','beg','end','gene_name','gene_id','strand','gene_type']
			# for genes with duplicated records, only keep the longest gene
			df_gene['length']=df_gene.end - df_gene.beg
			df_gene.sort_values('length',ascending=False,inplace=True) # type: ignore
			df_gene.drop_duplicates('gene_symbol',keep='first',inplace=True) # type: ignore
			df_gene.set_index('gene_symbol',inplace=True)
			for col in ['chrom','beg','end','strand','gene_type','gene_id','length']:
				adata.var[col]=adata.var_names.map(df_gene[col].to_dict())

		if normalization=='CPM':
			# for new sc-RNA-seq pipeline, CPM is equal to TPM?
			sc.pp.normalize_total(adata, target_sum=target_sum)
			sc.pp.log1p(adata) # log(CPM)
			adata.uns['Normalization']='log(CPM)'
		else: #TPM
			assert not gtf is None, "For TPM normalization, please provide gtf file."
			adata=cal_tpm(adata,target_sum=target_sum,length_fillna=1000)
	if not save is None:
		outdir=os.path.dirname(os.path.abspath(os.path.expanduser(save)))
		if not os.path.exists(outdir):
			os.makedirs(outdir,exist_ok=True)
		outfile=os.path.expanduser(save)
		adata.write_h5ad(outfile)
	else:
		return adata
