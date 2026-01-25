import os
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from .utils import normalize_mc_by_cell,parse_gtf
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

def cal_stats(adata_path,obs1,modality="RNA",expression_cutoff=0,
			  use_raw=False,normalize_per_cell=True,
			  clip_norm_value=10,sum_only=False):
	raw_adata=anndata.read_h5ad(os.path.expanduser(adata_path),backed='r')
	adata=raw_adata[obs1.index.tolist(),:].to_memory()
	raw_adata.file.close()
	if modality!='RNA':
		adata = normalize_mc_by_cell(
			use_adata=adata, normalize_per_cell=normalize_per_cell,
			clip_norm_value=clip_norm_value,
			hypo_score=False,verbose=0)
	else:
		if use_raw and not adata.raw is None:
			# adata.X=adata.raw.X.copy()
			adata_raw=adata.raw[:,adata.var_names.tolist()].to_adata()
			adata.X=adata_raw[adata.obs_names.tolist(),adata.var_names.tolist()].X.copy() # type: ignore
			del adata_raw
	df_data=adata.to_df() # rows are cells, columns are genes
	# Compute per-gene statistics (min, q25, q50, q75, max, sum) across cells.
	# Use NumPy's nanpercentile and nansum which are fast (uses quickselect under the hood).
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
	if sum_only:
		return sums,frac,df_data.columns.tolist()
	qs = np.nanpercentile(df_data.values, [0, 25, 50, 75, 100], axis=0)
	std=np.nanstd(df_data.values,axis=0)
	return qs,sums,std,frac,df_data.columns.tolist()

def cal_tpm(adata,target_sum=1e6,length_fillna=1000):
	assert 'length' in adata.var.columns.tolist(), "For TPM normalization, gene length information is required in adata.var['length'], please provide gene_meta"
	adata.var.length.fillna(length_fillna,inplace=True)
	# RPK (Reads Per Kilobase)
	counts = adata.to_df() #row are cell types and columns are genes
	# if hasattr(counts, 'toarray'):
	# 	counts = counts.toarray()  # Convert sparse to dense if needed
	lengths_kb = (adata.var['length'] / 1000).apply(lambda x:max(x,1)).to_dict() # keys are genes
	# RPK: divide each gene's counts by its length in kb
	rpk=counts.apply(lambda x: x / lengths_kb.get(x.name, 1)) # per column (gene), dataframe: rows are cell types, columns are genes
	# Calculate the "Per Million" Scaling Factor, Per-cell scaling factor: sum of RPK per cell
	rpk_sum = rpk.sum(axis=1).to_dict()  # Sum RPKs per cell (row); keys are cell types
	# TPM = (RPK / per_cell_sum) * 1e6
	tpm=rpk.apply(lambda x: (x / rpk_sum.get(x.name, 1)) * target_sum, axis=1) # Scale to TPM, for each row (cell type)
	# Store TPM in adata.layers
	adata.X = tpm.apply(np.log1p).values  # log(TPM)
	adata.uns['Normalization']='log(TPM)'
	return adata

def scrna2pseudobulk(
	adata_path,downsample=2000,
	obs_path=None,groupby="Group",use_raw=True,
	n_jobs=1,normalization=None,target_sum=1e6,gtf=None,save=None
):
	assert use_raw == True, "For normalization (CPM or TPM), please set use_raw=True"
	# assert modality=='RNA': # methylation
	raw_adata=anndata.read_h5ad(os.path.expanduser(adata_path),backed='r')
	if not obs_path is None:
		if isinstance(obs_path,str):
			obs=pd.read_csv(os.path.expanduser(obs_path),
				sep='\t',index_col=0)
		else:
			obs=obs_path.copy()
		overlapped_cells=list(set(raw_adata.obs_names.tolist()) & set(obs.index.tolist()))
		obs=obs.loc[overlapped_cells]
	else:
		obs=raw_adata.obs.copy()
		# raw_adata.obs[groupby]=raw_adata.obs.index.to_series().map(obs[groupby].to_dict())
	raw_adata.file.close()
	obs=obs.loc[obs[groupby].notna()]
	if not downsample is None:
		all_cells = obs.groupby(groupby).apply(
					lambda x: x.sample(downsample).index.tolist() if x.shape[0] > downsample else x.index.tolist()).sum()
	else:
		all_cells=obs.index.tolist()
	obs=obs.loc[all_cells]
	data={}
	if n_jobs==-1:
		n_jobs=os.cpu_count()
	with ProcessPoolExecutor(n_jobs) as executor:
		futures = {}
		for group in obs[groupby].unique():
			obs1=obs.loc[obs[groupby]==group]
			if obs1.shape[0]==0:
				continue
			future = executor.submit(
				cal_stats,adata_path=adata_path,obs1=obs1,
				use_raw=use_raw,sum_only=True
			)
			futures[future] = group
		logger.debug(f"Submitted {len(futures)} groups for pseudobulk calculation.")
		for future in as_completed(futures):
			group = futures[future]
			logger.debug(group)
			sums,frac,header = future.result()
			if 'sum' not in data:
				data['sum'] = []
			data['sum'].append(pd.Series(sums, name=group, index=header))
			frac.name = group
			if 'frac' not in data:
				data['frac'] = []
			data['frac'].append(pd.Series(frac, name=group,index=header))
	raw_adata.file.close()
	X=pd.concat(data['sum'],axis=1).T # sum of raw counts or normalized methylation fraction
	vc=raw_adata.obs.loc[all_cells][groupby].value_counts().to_frame(name='cell_count')
	adata = anndata.AnnData(X=X,obs=vc.loc[X.index.tolist()]) # put sum into adata.X
	adata.layers['frac']=pd.concat(objs=data['frac'],axis=1).T
	del data
	adata.raw=adata.copy()

	if not normalization is None and use_raw:
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
	vc_dict=vc.to_dict()['cell_count']
	adata.layers['mean']=adata.to_df().apply(lambda x:x/vc_dict[x.name],axis=1)
	if not save is None:
		outdir=os.path.dirname(os.path.abspath(os.path.expanduser(save)))
		if not os.path.exists(outdir):
			os.makedirs(outdir,exist_ok=True)
		outfile=os.path.expanduser(save)
		adata.write_h5ad(outfile)
	else:
		return adata

def stat_pseudobulk(
	adata_path,downsample=2000,
	obs_path=None,groupby="Group",use_raw=False,expression_cutoff=0,
	modality="RNA",n_jobs=1,normalize_per_cell=True,clip_norm_value=10,
	save=None
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
	else:
		obs=raw_adata.obs.copy()
		# raw_adata.obs[groupby]=raw_adata.obs.index.to_series().map(obs[groupby].to_dict())
	raw_adata.file.close()
	obs=obs.loc[obs[groupby].notna()]
	if not downsample is None:
		all_cells = obs.groupby(groupby).apply(
					lambda x: x.sample(downsample).index.tolist() if x.shape[0] > downsample else x.index.tolist()).sum()
	else:
		all_cells=obs.index.tolist()
	obs=obs.loc[all_cells]
	data={}
	if n_jobs==-1:
		n_jobs=os.cpu_count()
	with ProcessPoolExecutor(n_jobs) as executor:
		futures = {}
		for group in obs[groupby].unique():
			obs1=obs.loc[obs[groupby]==group]
			if obs1.shape[0]==0:
				continue
			future = executor.submit(
				cal_stats,adata_path=adata_path,obs1=obs1,modality=modality,
				expression_cutoff=expression_cutoff,
				use_raw=use_raw,normalize_per_cell=normalize_per_cell,
				clip_norm_value=clip_norm_value
			)
			futures[future] = group
		logger.debug(f"Submitted {len(futures)} groups for pseudobulk calculation.")
		for future in as_completed(futures):
			group = futures[future]
			logger.debug(group)
			qs,sums,std,frac,header = future.result()
			for k,v in zip(['min', 'q25', 'q50', 'q75', 'max', 'sum','std'], qs.tolist() + [sums.tolist(),std.tolist()]):
				if k not in data:
					data[k] = []
				data[k].append(pd.Series(v, name=group, index=header))
			frac.name = group
			if 'frac' not in data:
				data['frac'] = []
			data['frac'].append(pd.Series(frac, name=group,index=header))
	X=pd.concat(data['sum'],axis=1).T # sum of raw counts or normalized methylation fraction
	vc=obs[groupby].value_counts().to_frame(name='cell_count')
	adata = anndata.AnnData(X=X,obs=vc.loc[X.index.tolist()]) # put sum into adata.X
	for k in data:
		if k=='sum':
			continue
		adata.layers[k]=pd.concat(objs=data[k],axis=1).T
	del data
	vc_dict=vc.to_dict()['cell_count']
	adata.layers['mean']=adata.to_df().apply(lambda x:x/vc_dict[x.name],axis=1)
	if not save is None:
		outdir=os.path.dirname(os.path.abspath(os.path.expanduser(save)))
		if not os.path.exists(outdir):
			os.makedirs(outdir,exist_ok=True)
		outfile=os.path.expanduser(save)
		adata.write_h5ad(outfile)
	else:
		return adata

def export_pseudobulk_adata(adata,outdir,use_raw):
	"""
	Export pseudobulk adata to bed
	"""
	outdir=os.path.expanduser(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir,exist_ok=True)
	if not os.path.exists(outdir):
		os.makedirs(outdir,exist_ok=True)
	if isinstance(adata,str):
		adata=anndata.read_h5ad(os.path.expanduser(adata))
	else:
		adata=adata
	if use_raw:
		data=adata.raw.to_adata().to_df().T # raw counts
	else:
		data=adata.to_df().T # CPM, log(CPM) or ..
	if "chrom" in adata.var.columns.tolist():
		data.insert(0,"chrom",adata.var.loc[data.index.tolist(),"chrom"].tolist())
	else:
		data.insert(0,"chrom",data.index.to_series().apply(lambda x:x.split(':')[0]))
	if "start" in adata.var.columns.tolist():
		data.insert(1,"start",adata.var.loc[data.index.tolist(),"start"].tolist())
	else:
		data.insert(1,"start",data.index.to_series().apply(lambda x:x.split(':')[1].split('-')[0]))
	if "end" in adata.var.columns.tolist():
		data.insert(2,"end",adata.var.loc[data.index.tolist(),"end"].tolist())
	else:
		data.insert(2,"end",data.index.to_series().apply(lambda x:x.split(':')[1].split('-')[1]))
	data.insert(3,"features",data.index.tolist())
	if "strand" in adata.var.columns.tolist():
		data.insert(4,"strand",adata.var.loc[data.index.tolist(),"strand"].tolist())
	else:
		data.insert(4,"strand","+")
	data=data.loc[(data.chrom.notna()) & (data.start.notna()) & (data.end.notna())]
	data.start=data.start.astype(int)
	data.end=data.end.astype(int)
	data.sort_values(['chrom','start','end'],ascending=True,inplace=True)
	for col in data.columns.tolist()[4:]:
		data.loc[:,['chrom','start','end','features',col,'strand']].to_csv(
			os.path.join(outdir,f"{col.replace(' ','_')}.bed"),
				sep='\t',index=False,header=False)
		df=data.loc[:,['features',col]]
		df.to_csv(os.path.join(outdir,f"{col.replace(' ','_')}.txt"),
			sep='\t',index=False,header=False)

def downsample_adata(adata_path,groupby="Group",obs_path=None,
					 outfile="Group.downsample_1500.h5ad",
					 downsample=1500):
	adata_path=os.path.expanduser(adata_path)
	outfile=os.path.expanduser(outfile)
	adata=anndata.read_h5ad(adata_path,backed='r')
	if not obs_path is None:
		if isinstance(obs_path,str):
			obs=pd.read_csv(os.path.expanduser(obs_path),
				sep='\t',index_col=0)
		else:
			obs=obs_path.copy()
		overlapped_cells=list(set(adata.obs_names.tolist()) & set(obs.index.tolist()))
		obs=obs.loc[overlapped_cells]
	else:
		obs=adata.obs.copy()
	keep_cells = obs.loc[obs[groupby].notna()].groupby(groupby).apply(
		lambda x: x.sample(downsample).index.tolist() if x.shape[0] > downsample else x.index.tolist()).sum()
	adata[keep_cells,:].write_h5ad(outfile,compression='gzip')
	adata.file.close()
