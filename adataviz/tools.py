import os
import numpy as np
import pandas as pd
import anndata
from .plotting import categorical_scatter,continuous_scatter

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

def parse_gtf(gtf="gencode.v43.annotation.gtf"):
    df=pd.read_csv(os.path.expanduser(gtf),sep='\t',header=None,
                    comment="#",usecols=[0,2,3,4,6,8],
                    names=['chrom','record_type','beg','end','strand','information'])
    cols=['gene_id','gene_type','gene_name']
    def parse_info(x):
        x=x.replace('"','')
        D={}
        for item in x.strip().rstrip(';').split(';'):
            k,v=item.strip().split(' ')
            D[k.strip()]=v.strip()
        return D

    df['info_dict']=df.information.apply(parse_info)
    for col in cols:
        df[col]=df.info_dict.apply(lambda x:x.get(col,''))
    df=df.loc[:,['chrom','beg','end','gene_name','gene_id','strand','gene_type']].drop_duplicates()
    return df # 'chrom','start','end','gene_symbol','strand','gene_type'

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
