import os
import numpy as np
import pandas as pd
import anndata
from .utils import normalize_mc_by_cell, parse_gtf
from concurrent.futures import ProcessPoolExecutor, as_completed
from loguru import logger as logger


def merge_adata_regions(
    pseudobulk_adata_path,
    bin_size=5000,
    use_raw=True,
    res=100000,
    filter_chroms=True,
    boundary=None,
    exclude_chroms=["chrY", "chrM", "chrX"],
):
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
        tmp_name = ""
        for x, idx, name in zip(list_x, indexes, names):
            chrom = name.split("_")[0]
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

    import scanpy as sc

    adata = load_adata(pseudobulk_adata_path, backed=None)
    if use_raw:
        if adata.raw is not None:
            adata = adata.raw.to_adata()
        else:
            logger.warning("adata.raw is None!!")
    data = adata.to_df().T
    fea = "100kb" if res == 100000 else "10kb" if res == 10000 else "25kb"
    if boundary is None:
        boundary = True if fea == "25kb" else False
    groups = data.columns.tolist()
    _idx = data.index.to_series()
    data["chrom"] = _idx.str.split(":").str[0]
    data["start"] = _idx.str.split(":").str[1].str.split("-").str[0].astype(int)
    data["BinID"] = (data["start"] // bin_size).astype(int)

    if not boundary:
        # merge 5kb into 100kb or 10kb
        data[fea] = data["chrom"] + "_" + (data["start"] // res).astype(str)
        # index, for example, chr1_0, chr1_0, chr1_0...,chr1_1
        data = data.loc[:, groups + [fea]].groupby(fea).sum().T
    else:  # for 25kb, get the domain doundary (10kb), not the 25kb windows.
        logger.info("Generating results for boundaries of 25kb bins")
        ids = data["BinID"] % 5  # [0,1],2,[3,4,0,1],2,[3,4,0,1],...
        data[fea] = data["chrom"] + "_" + (data["BinID"] // 5).astype(str)
        name2index = assign_boundary(
            data.index.tolist(), data[fea].tolist(), ids.tolist()
        )
        idx2name = {}
        for name in name2index:
            for idx in name2index[name]:
                idx2name[idx] = name
        data["NAME"] = data.index.to_series().map(idx2name)
        data = data.loc[data["NAME"].notna()]
        data = data.loc[data["NAME"].apply(lambda x: x.split("_")[0]) == data.chrom]
        data = data.loc[:, groups + ["NAME"]].groupby("NAME").sum().T
    if use_raw:  # convert sum of raw counts into CPM
        adata = anndata.AnnData(X=data)
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)  # log(CPM)
        data = adata.to_df()
    if filter_chroms:
        s_col = data.columns.to_frame()
        chrom_s = s_col.index.to_series().str.split("_").str[0]
        s_col["chrom"] = chrom_s
        keep_cols = s_col.loc[
            ~chrom_s.isin(exclude_chroms) & (chrom_s.str.len() < 6)
        ].index.tolist()
        # use_rows = list(set(mc_df.index.tolist()) & set(cellclass2majortype[group]))
        data = data.loc[:, keep_cols]
    # rows are 100kb bins, columns are cell types
    return data  # to do next: subset rows (df_bin.index) and columns (cell types order)


def cal_stats(
    adata_path,
    obs1,
    modality="RNA",
    expression_cutoff=0,
    use_raw=False,
    normalize_per_cell=True,
    clip_norm_value=10,
    sum_only=False,
):
    """Compute per-gene summary statistics for a subset of cells.

    Loads the cells identified by *obs1* from a backed ``.h5ad`` file and
    computes column-wise (per-gene) statistics including percentiles, sum,
    standard deviation, and the fraction of cells expressing each gene.

    Parameters
    ----------
    adata_path : str
            Path to the ``.h5ad`` file.
    obs1 : pandas.DataFrame
            Subset of cell metadata whose index selects the cells to use.
    modality : str, optional
            Data modality (``'RNA'``, ``'ATAC'``, or methylation).
            Default is ``'RNA'``.
    expression_cutoff : float, optional
            Threshold above which a gene is considered expressed (for RNA/ATAC).
            Default is ``0``.
    use_raw : bool, optional
            Whether to use ``adata.raw`` for the expression matrix.
            Default is ``False``.
    normalize_per_cell : bool, optional
            Whether to normalise methylation per cell. Default is ``True``.
    clip_norm_value : float, optional
            Clipping value for per-cell normalisation. Default is ``10``.
    sum_only : bool, optional
            If ``True``, return only sums and fractions (faster for pseudobulk
            creation). Default is ``False``.

    Returns
    -------
    tuple
            If ``sum_only=True``: ``(sums, frac, gene_names)``.
            Otherwise: ``(qs, sums, std, frac, gene_names)`` where *qs* is a
            ``(5, n_genes)`` array of ``[min, q25, q50, q75, max]``.
    """
    raw_adata = load_adata(adata_path, backed="r")
    adata = raw_adata[obs1.index.tolist(), :].to_memory()
    raw_adata.file.close()
    if modality not in ["RNA", "ATAC"]:
        adata = normalize_mc_by_cell(
            use_adata=adata,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
            hypo_score=False,
            verbose=0,
        )
    else:
        if use_raw and adata.raw is not None:
            # adata.X=adata.raw.X.copy()
            adata_raw = adata.raw[:, adata.var_names.tolist()].to_adata()
            adata.X = adata_raw[
                adata.obs_names.tolist(), adata.var_names.tolist()
            ].X.copy()  # type: ignore
            del adata_raw
    df_data = adata.to_df()  # rows are cells, columns are genes
    # Compute per-gene statistics (min, q25, q50, q75, max, sum) across cells.
    # Use NumPy's nanpercentile and nansum which are fast (uses quickselect under the hood).
    sums = np.nansum(df_data.values, axis=0)  # for each column
    # fraction of cells expressing (or hypomethylated) the gene
    if modality not in ["RNA", "ATAC"]:  # methylation, cutoff = 1
        # frac = df_data.apply(lambda x: x[x < 1].shape[0] / x.shape[0])
        # vectorized: count values < 1 per column divided by number of cells
        frac = (df_data < 1).sum(axis=0) / float(df_data.shape[0])
    else:  # for RNA
        # frac = df_data.apply(lambda x: x[x > expression_cutoff].shape[0] / x.shape[0])
        # vectorized: count values > cutoff per column divided by number of cells
        frac = (df_data > expression_cutoff).sum(axis=0) / float(df_data.shape[0])
    if sum_only:
        return sums, frac, df_data.columns.tolist()
    qs = np.nanpercentile(df_data.values, [0, 25, 50, 75, 100], axis=0)
    std = np.nanstd(df_data.values, axis=0)
    return qs, sums, std, frac, df_data.columns.tolist()


def cal_tpm(adata, target_sum=1e6, length_fillna=1000):
    """Normalise an AnnData object to Transcripts Per Million (TPM).

    Requires ``adata.var['length']`` to contain gene lengths. The result
    is stored as ``log1p(TPM)`` in ``adata.X``.

    Parameters
    ----------
    adata : AnnData
            Annotated data matrix with raw counts in ``.X`` and a ``'length'``
            column in ``.var``.
    target_sum : float, optional
            Target total per cell after scaling. Default is ``1e6``.
    length_fillna : int, optional
            Fill value for missing gene lengths. Default is ``1000``.

    Returns
    -------
    AnnData
            The same object with ``adata.X`` replaced by ``log1p(TPM)``.
    """
    assert "length" in adata.var.columns.tolist(), (
        "For TPM normalization, gene length information is required in adata.var['length'], please provide gene_meta"
    )
    adata.var.length.fillna(length_fillna, inplace=True)
    # RPK (Reads Per Kilobase)
    counts = adata.to_df()  # row are cell types and columns are genes
    lengths_kb = (adata.var["length"] / 1000).clip(lower=1)  # Series, index=genes
    # RPK: divide each gene's counts by its length in kb (vectorized broadcasting)
    rpk = counts.div(lengths_kb, axis=1)
    # Calculate the "Per Million" Scaling Factor, Per-cell scaling factor: sum of RPK per cell
    rpk_sum = rpk.sum(axis=1)  # Series, index=cell types
    # TPM = (RPK / per_cell_sum) * 1e6 (vectorized broadcasting)
    tpm = rpk.div(rpk_sum, axis=0) * target_sum
    # Store TPM in adata.layers
    adata.X = np.log1p(tpm.values)  # log(TPM)
    adata.uns["Normalization"] = "log(TPM)"
    return adata


def scrna2pseudobulk(
    adata_path,
    downsample=2000,
    obs_path=None,
    groupby="Group",
    use_raw=True,
    n_jobs=1,
    normalization=None,
    target_sum=1e6,
    gtf=None,
    save=None,
):
    """Aggregate single-cell RNA data into pseudobulk profiles per group.

    For each unique value in *groupby*, sums raw counts across cells
    (optionally downsampled), then optionally normalises to CPM or TPM.

    Parameters
    ----------
    adata_path : str
            Path to the ``.h5ad`` file.
    downsample : int or None, optional
            Maximum number of cells per group to use. ``None`` disables
            downsampling. Default is ``2000``.
    obs_path : str, DataFrame, or None, optional
            Cell metadata. Path to a TSV file or a DataFrame. If ``None``,
            uses ``adata.obs``.
    groupby : str, optional
            Column defining groups. Default is ``'Group'``.
    use_raw : bool, optional
            Whether to use raw counts. Must be ``True`` for normalisation.
            Default is ``True``.
    n_jobs : int, optional
            Number of parallel workers. ``-1`` uses all CPUs. Default is ``1``.
    normalization : str or None, optional
            Normalisation method: ``'CPM'`` or ``'TPM'``. ``None`` keeps raw
            sums. Default is ``None``.
    target_sum : float, optional
            Target total for CPM/TPM scaling. Default is ``1e6``.
    gtf : str or None, optional
            Path to a GTF file (required for TPM normalisation).
    save : str or None, optional
            Path to save the resulting ``.h5ad``. If ``None``, the AnnData
            object is returned.

    Returns
    -------
    AnnData or None
            Pseudobulk AnnData if ``save`` is ``None``.
    """
    assert use_raw, "For normalization (CPM or TPM), please set use_raw=True"
    # assert modality=='RNA': # methylation
    raw_adata = load_adata(adata_path, backed="r")
    if obs_path is not None:
        if isinstance(obs_path, str):
            obs = pd.read_csv(os.path.expanduser(obs_path), sep="\t", index_col=0)
        else:
            obs = obs_path.copy()
        overlapped_cells = list(
            set(raw_adata.obs_names.tolist()) & set(obs.index.tolist())
        )
        obs = obs.loc[overlapped_cells]
    else:
        obs = raw_adata.obs.copy()
        # raw_adata.obs[groupby]=raw_adata.obs.index.to_series().map(obs[groupby].to_dict())
    raw_adata.file.close()
    obs = obs.loc[obs[groupby].notna()]
    if downsample is not None:
        all_cells = (
            obs.groupby(groupby)
            .apply(
                lambda x: x.sample(downsample).index.tolist()
                if x.shape[0] > downsample
                else x.index.tolist()
            )
            .sum()
        )
    else:
        all_cells = obs.index.tolist()
    obs = obs.loc[all_cells]
    data = {}
    if n_jobs == -1:
        n_jobs = os.cpu_count()
    with ProcessPoolExecutor(n_jobs) as executor:
        futures = {}
        for group in obs[groupby].unique():
            obs1 = obs.loc[obs[groupby] == group]
            if obs1.shape[0] == 0:
                continue
            future = executor.submit(
                cal_stats,
                adata_path=adata_path,
                obs1=obs1,
                use_raw=use_raw,
                sum_only=True,
            )
            futures[future] = group
        logger.debug(f"Submitted {len(futures)} groups for pseudobulk calculation.")
        for future in as_completed(futures):
            group = futures[future]
            logger.debug(group)
            sums, frac, header = future.result()
            if "sum" not in data:
                data["sum"] = []
            data["sum"].append(pd.Series(sums, name=group, index=header))
            frac.name = group
            if "frac" not in data:
                data["frac"] = []
            data["frac"].append(pd.Series(frac, name=group, index=header))
    raw_adata.file.close()
    X = pd.concat(
        data["sum"], axis=1
    ).T  # sum of raw counts or normalized methylation fraction
    vc = (
        raw_adata.obs.loc[all_cells][groupby].value_counts().to_frame(name="cell_count")
    )
    adata = anndata.AnnData(X=X, obs=vc.loc[X.index.tolist()])  # put sum into adata.X
    adata.layers["frac"] = pd.concat(objs=data["frac"], axis=1).T
    del data
    adata.raw = adata.copy()

    if normalization is not None and use_raw:
        # Calculate CPM or TPM only if aggfunc is sum
        logger.info(f"Normalizing pseudobulk adata using {normalization} method.")
        if gtf is not None:
            df_gene = parse_gtf(gtf=gtf)
            # ['chrom','beg','end','gene_name','gene_id','strand','gene_type']
            # for genes with duplicated records, only keep the longest gene
            df_gene["length"] = df_gene.end - df_gene.beg
            df_gene.sort_values("length", ascending=False, inplace=True)  # type: ignore
            df_gene.drop_duplicates("gene_name", keep="first", inplace=True)  # type: ignore
            df_gene.set_index("gene_name", inplace=True)
            for col in [
                "chrom",
                "beg",
                "end",
                "strand",
                "gene_type",
                "gene_id",
                "length",
            ]:
                adata.var[col] = adata.var_names.map(df_gene[col].to_dict())

        if normalization == "CPM":
            import scanpy as sc

            # for new sc-RNA-seq pipeline, CPM is equal to TPM?
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)  # log(CPM)
            adata.uns["Normalization"] = "log(CPM)"
        else:  # TPM
            assert gtf is not None, "For TPM normalization, please provide gtf file."
            adata = cal_tpm(adata, target_sum=target_sum, length_fillna=1000)
    cell_counts = vc["cell_count"].reindex(adata.obs_names)
    adata.layers["mean"] = adata.to_df().div(cell_counts, axis=0)
    if save is not None:
        outdir = os.path.dirname(os.path.abspath(os.path.expanduser(save)))
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        outfile = os.path.expanduser(save)
        adata.write_h5ad(outfile)
    else:
        return adata


def stat_pseudobulk(
    adata_path,
    downsample=2000,
    obs_path=None,
    groupby="Group",
    use_raw=False,
    expression_cutoff=0,
    modality="RNA",
    n_jobs=1,
    normalize_per_cell=True,
    clip_norm_value=10,
    save=None,
):
    """Compute pseudobulk summary statistics (min/q25/q50/q75/max/std/frac) per group.

    For each group defined by *groupby*, computes per-gene percentiles,
    standard deviation, and the fraction of expressing (or hypomethylated)
    cells, storing results in ``adata.layers``.

    Parameters
    ----------
    adata_path : str
            Path to the ``.h5ad`` file.
    downsample : int or None, optional
            Maximum number of cells per group. Default is ``2000``.
    obs_path : str, DataFrame, or None, optional
            Cell metadata. If ``None``, uses ``adata.obs``.
    groupby : str, optional
            Grouping column. Default is ``'Group'``.
    use_raw : bool, optional
            Whether to use ``adata.raw``. Default is ``False``.
    expression_cutoff : float, optional
            Threshold for marking a gene as expressed. Default is ``0``.
    modality : str, optional
            Data modality (``'RNA'``, ``'ATAC'``, or methylation).
            Default is ``'RNA'``.
    n_jobs : int, optional
            Number of parallel workers. ``-1`` uses all CPUs. Default is ``1``.
    normalize_per_cell : bool, optional
            Whether to normalise methylation per cell. Default is ``True``.
    clip_norm_value : float, optional
            Clipping value for normalisation. Default is ``10``.
    save : str or None, optional
            Path to save the resulting ``.h5ad``. If ``None``, the AnnData
            object is returned.

    Returns
    -------
    AnnData or None
            Pseudobulk AnnData with stat layers if ``save`` is ``None``.
    """
    if modality not in ["RNA", "ATAC"]:  # methylation
        assert normalize_per_cell, "For methylation, normalize_per_cell should be True"
    raw_adata = load_adata(adata_path, backed="r")
    if obs_path is not None:
        if isinstance(obs_path, str):
            obs = pd.read_csv(os.path.expanduser(obs_path), sep="\t", index_col=0)
        else:
            obs = obs_path.copy()
        overlapped_cells = list(
            set(raw_adata.obs_names.tolist()) & set(obs.index.tolist())
        )
        obs = obs.loc[overlapped_cells]
    else:
        obs = raw_adata.obs.copy()
        # raw_adata.obs[groupby]=raw_adata.obs.index.to_series().map(obs[groupby].to_dict())
    raw_adata.file.close()
    obs = obs.loc[obs[groupby].notna()]
    if downsample is not None:
        all_cells = (
            obs.groupby(groupby)
            .apply(
                lambda x: x.sample(downsample).index.tolist()
                if x.shape[0] > downsample
                else x.index.tolist()
            )
            .sum()
        )
    else:
        all_cells = obs.index.tolist()
    obs = obs.loc[all_cells]
    data = {}
    if n_jobs == -1:
        n_jobs = os.cpu_count()
    with ProcessPoolExecutor(n_jobs) as executor:
        futures = {}
        for group in obs[groupby].unique():
            obs1 = obs.loc[obs[groupby] == group]
            if obs1.shape[0] == 0:
                continue
            future = executor.submit(
                cal_stats,
                adata_path=adata_path,
                obs1=obs1,
                modality=modality,
                expression_cutoff=expression_cutoff,
                use_raw=use_raw,
                normalize_per_cell=normalize_per_cell,
                clip_norm_value=clip_norm_value,
            )
            futures[future] = group
        logger.debug(f"Submitted {len(futures)} groups for pseudobulk calculation.")
        for future in as_completed(futures):
            group = futures[future]
            logger.debug(group)
            qs, sums, std, frac, header = future.result()
            for k, v in zip(
                ["min", "q25", "q50", "q75", "max", "sum", "std"],
                qs.tolist() + [sums.tolist(), std.tolist()],
            ):
                if k not in data:
                    data[k] = []
                data[k].append(pd.Series(v, name=group, index=header))
            frac.name = group
            if "frac" not in data:
                data["frac"] = []
            data["frac"].append(pd.Series(frac, name=group, index=header))
    X = pd.concat(
        data["sum"], axis=1
    ).T  # sum of raw counts or normalized methylation fraction
    vc = obs[groupby].value_counts().to_frame(name="cell_count")
    adata = anndata.AnnData(X=X, obs=vc.loc[X.index.tolist()])  # put sum into adata.X
    for k in data:
        if k == "sum":
            continue
        adata.layers[k] = pd.concat(objs=data[k], axis=1).T
    del data
    cell_counts = vc["cell_count"].reindex(adata.obs_names)
    adata.layers["mean"] = adata.to_df().div(cell_counts, axis=0)
    if save is not None:
        outdir = os.path.dirname(os.path.abspath(os.path.expanduser(save)))
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        outfile = os.path.expanduser(save)
        adata.write_h5ad(outfile)
    else:
        return adata


def normalize_adata(
    adata,
    embedding=True,
    outfile=None,
    n_top_features=5000,
    n_pcs=50,
    min_cells=5,
    batch_col=None,
    normalization="CPM",
    gtf=None,
    flanking=None,
    target_sum=1e4,
):
    """Normalise an AnnData and optionally compute embeddings.

    Performs gene filtering, CPM or TPM normalisation, and optionally
    highly-variable-gene selection, PCA, Harmony batch correction,
    neighbour graph, Leiden clustering, and UMAP.

    Parameters
    ----------
    adata : AnnData or str
            Annotated data matrix or path to an ``.h5ad`` file.
    embedding : bool, optional
            Whether to compute PCA, neighbours, Leiden, and UMAP.
            Default is ``True``.
    outfile : str or None, optional
            Path to save the processed ``.h5ad``. If ``None``, the AnnData
            object is returned.
    n_top_features : int, optional
            Number of highly variable genes to select. Default is ``5000``.
    n_pcs : int, optional
            Number of principal components. Default is ``50``.
    min_cells : int, optional
            Minimum number of cells a gene must be expressed in.
            Default is ``5``.
    batch_col : str or None, optional
            Column in ``adata.obs`` for Harmony batch correction.
            If ``None``, no batch correction is performed.
    normalization : str, optional
            Normalisation method: ``'CPM'`` or ``'TPM'``. Default is
            ``'CPM'``.
    gtf : str or None, optional
            Path to a GTF file (required for TPM).
    flanking : int or None, optional
            Base pairs added to gene body for TPM length calculation.
    target_sum : float, optional
            Target total per cell for scaling. Default is ``1e4``.

    Returns
    -------
    AnnData or None
            Normalised AnnData if ``outfile`` is ``None``.
    """
    import scanpy as sc

    adata = load_adata(adata, backed=None)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    if normalization == "CPM":
        # for new sc-RNA-seq pipeline, CPM is equal to TPM?
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)  # log(CPM)
        adata.uns["Normalization"] = "log(CPM)"
    else:  # TPM
        assert gtf is not None, "For TPM normalization, please provide gtf file."
        df_gene = parse_gtf(gtf=gtf)
        # ['chrom','beg','end','gene_name','gene_id','strand','gene_type']
        # for genes with duplicated records, only keep the longest gene
        if flanking is not None:
            logger.info(
                f"Adding flanking {flanking} bp to gene body for TPM calculation."
            )
            df_gene["beg"] = df_gene.beg.apply(lambda x: max(x - flanking, 0))
            df_gene["end"] = df_gene.end + flanking
        df_gene["length"] = df_gene.end - df_gene.beg
        df_gene.sort_values("length", ascending=False, inplace=True)  # type: ignore
        df_gene.drop_duplicates("gene_name", keep="first", inplace=True)  # type: ignore
        df_gene.set_index("gene_name", inplace=True)
        for col in ["chrom", "beg", "end", "strand", "gene_type", "gene_id", "length"]:
            adata.var[col] = adata.var_names.map(df_gene[col].to_dict())
        adata = cal_tpm(adata, target_sum=target_sum, length_fillna=1000)

    if embedding:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_features)
        sc.tl.pca(adata)
        sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs)  # log=True
        if batch_col is not None:
            import scanpy.external as sce

            sce.pp.harmony_integrate(
                adata, key=batch_col, basis="X_pca", max_iter_harmony=50
            )
            use_rep = "X_pca_harmony"
        else:
            use_rep = "X_pca"  #'X_pca_harmony'
        sc.pp.neighbors(adata, use_rep=use_rep)
        sc.tl.leiden(adata, resolution=1)
        sc.tl.umap(adata)
    if outfile is not None:
        outfile = os.path.expanduser(outfile)
        outdir = os.path.dirname(os.path.abspath(outfile))
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        adata.write_h5ad(outfile)
    else:
        return adata


def export_pseudobulk_adata(adata, outdir="pseudobulk.bed", use_raw=False):
    """
    Export pseudobulk adata to bed
    """
    outdir = os.path.expanduser(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    adata = load_adata(adata, backed=None)
    if use_raw:
        data = adata.raw.to_adata().to_df().T  # raw counts
    else:
        data = adata.to_df().T  # CPM, log(CPM) or ..
    if "chrom" in adata.var.columns.tolist():
        data.insert(0, "chrom", adata.var.loc[data.index.tolist(), "chrom"].tolist())
    else:
        data.insert(0, "chrom", data.index.to_series().apply(lambda x: x.split(":")[0]))
    beg = "beg" if "beg" in adata.var.columns.tolist() else "start"
    if beg in adata.var.columns.tolist():
        data.insert(1, "start", adata.var.loc[data.index.tolist(), beg].tolist())
    else:
        data.insert(
            1,
            "start",
            data.index.to_series().apply(lambda x: x.split(":")[1].split("-")[0]),
        )
    if "end" in adata.var.columns.tolist():
        data.insert(2, "end", adata.var.loc[data.index.tolist(), "end"].tolist())
    else:
        data.insert(
            2,
            "end",
            data.index.to_series().apply(lambda x: x.split(":")[1].split("-")[1]),
        )
    data.insert(3, "features", data.index.tolist())
    if "strand" in adata.var.columns.tolist():
        data.insert(4, "strand", adata.var.loc[data.index.tolist(), "strand"].tolist())
    else:
        data.insert(4, "strand", "+")
    data = data.loc[(data.chrom.notna()) & (data.start.notna()) & (data.end.notna())]
    data.start = data.start.astype(int)
    data.end = data.end.astype(int)
    data.sort_values(["chrom", "start", "end"], ascending=True, inplace=True)
    for col in data.columns.tolist()[4:]:
        if col in ["chrom", "start", "beg", "end", "strand", "features"]:
            continue
        data.loc[:, ["chrom", "start", "end", "features", col, "strand"]].to_csv(
            os.path.join(outdir, f"{col.replace(' ', '_')}.bed"),
            sep="\t",
            index=False,
            header=False,
        )
        df = data.loc[:, ["features", col]]
        df.to_csv(
            os.path.join(outdir, f"{col.replace(' ', '_')}.txt"),
            sep="\t",
            index=False,
            header=False,
        )


def load_adata(adata, backed="r"):
    """Load an AnnData object from a file path or pass through an existing one.

    Parameters
    ----------
    adata : AnnData or str
            An existing AnnData object or a path to an ``.h5ad`` file.
    backed : str or None, optional
            Backing mode for ``anndata.read_h5ad()``. Default is ``'r'``
            (read-only backed mode).

    Returns
    -------
    AnnData
    """
    if isinstance(adata, str):
        adata = anndata.read_h5ad(os.path.expanduser(adata), backed=backed)
    else:
        assert isinstance(adata, anndata.AnnData)
    return adata


def load_obs(obs):
    """Load cell metadata from various sources.

    Accepts a file path (TSV, CSV, or ``.h5ad``) or a DataFrame. For
    ``.h5ad`` files, reads ``adata.obs`` in backed mode and closes the
    file handle.

    Parameters
    ----------
    obs : str or DataFrame
            Path to a metadata file or an existing DataFrame.

    Returns
    -------
    pandas.DataFrame
            A copy of the cell metadata.
    """
    if isinstance(obs, str) and not obs.endswith(".h5ad"):
        obs_path = os.path.abspath(os.path.expanduser(obs))
        sep = "\t" if obs_path.endswith(".tsv") or obs_path.endswith(".txt") else ","
        obs = pd.read_csv(obs_path, index_col=0, sep=sep)
        return obs
    elif isinstance(obs, str) and obs.endswith(".h5ad"):
        adata = anndata.read_h5ad(os.path.expanduser(obs), backed="r")
        obs = adata.obs.copy()
        adata.file.close()
        return obs
    else:
        assert isinstance(obs, pd.DataFrame)
        return obs.copy()


def load_color_palette(palette=None, adata=None, groups=[]):
    """Load a colour palette from an Excel file or from ``adata.uns``.

    Parameters
    ----------
    palette : str or None, optional
            Path to an Excel file where each sheet is named after a grouping
            column and contains a ``'Hex'`` column with colour codes.
            If ``None``, colours are extracted from ``adata.uns``.
    adata : AnnData or None, optional
            Annotated data matrix used to look up colours in ``adata.uns``
            when *palette* is not provided.
    groups : str or list of str, optional
            Grouping column name(s) to extract colours for.

    Returns
    -------
    dict or None
            If a single group is given, returns ``{category: hex_color}``.
            If multiple groups, returns ``{group_col: {category: hex_color}}``.
            Returns ``None`` if no colours are found.
    """
    # read color palette
    if isinstance(groups, str):
        groups = [groups]
    assert isinstance(groups, list)
    color_palette = {}
    if palette is not None and os.path.exists(os.path.expanduser(palette)):
        palette = os.path.abspath(os.path.expanduser(palette))
        D = pd.read_excel(palette, sheet_name=None, index_col=0)
        for group_col in groups:
            if group_col in D:
                color_palette[group_col] = D[group_col].Hex.to_dict()
            else:
                assert "+" in group_col, f"{group_col} not found in the palette file."
                for group in group_col.split("+"):
                    assert group in D, f"{group} not found in the palette file."
                    color_palette[group] = D[group].Hex.to_dict()
    else:
        if adata is None:
            return None
        # assert not adata is None, "If palette not provided, adata is required to extract color information from adata.obs"
        for group_col in groups:
            if f"{group_col}_colors" in adata.uns:
                group_series = adata.obs[group_col]
                if not pd.api.types.is_categorical_dtype(group_series):
                    group_series = group_series.astype("category")
                color_palette[group_col] = {
                    cluster: color
                    for cluster, color in zip(
                        group_series.cat.categories.tolist(),
                        adata.uns[f"{group_col}_colors"],
                    )
                }
    if len(color_palette) == 0:
        return None
    if len(groups) == 1:
        return color_palette[groups[0]]
    return color_palette


def get_obs(adata_path, add_coord=True, usecols=None, index_name="cell", outfile=None):
    """Extract ``adata.obs`` and optionally append embedding coordinates.

    Parameters
    ----------
    adata_path : str
            Path to an ``.h5ad`` file.
    add_coord : bool, optional
            Whether to append UMAP and/or t-SNE coordinates from
            ``adata.obsm``. Default is ``True``.
    usecols : list of str or None, optional
            Subset of columns to keep. ``None`` keeps all.
    index_name : str, optional
            Name assigned to the index column. Default is ``'cell'``.
    outfile : str or None, optional
            Path to save the result as a TSV file. If ``None``, the
            DataFrame is returned.

    Returns
    -------
    pandas.DataFrame or str
            The obs DataFrame (reset index) or the *outfile* path.
    """
    adata = load_adata(adata_path, backed="r")
    obs = adata.obs.copy()
    obs.index.name = index_name
    if add_coord:
        for coord in ["umap", "tsne"]:
            if f"X_{coord}" not in adata.obsm:
                continue
            df_coord = pd.DataFrame(
                adata.obsm[f"X_{coord}"],
                columns=[f"{coord}_0", f"{coord}_1"],
                index=adata.obs_names,
            )
            obs[f"{coord}_0"] = obs.index.to_series().map(
                df_coord[f"{coord}_0"].to_dict()
            )
            obs[f"{coord}_1"] = obs.index.to_series().map(
                df_coord[f"{coord}_1"].to_dict()
            )
    adata.file.close()
    if usecols is not None:
        obs = obs.loc[:, usecols]
    if outfile is not None:
        obs.to_csv(os.path.expanduser(outfile), sep="\t")
        return outfile
    return obs.reset_index()


def downsample_adata(
    adata_path,
    groupby="Group",
    obs_path=None,
    outfile="Group.downsample_1500.h5ad",
    downsample=1500,
):
    """Downsample an AnnData to at most *downsample* cells per group.

    Parameters
    ----------
    adata_path : str
            Path to the source ``.h5ad`` file.
    groupby : str, optional
            Column defining groups. Default is ``'Group'``.
    obs_path : str, DataFrame, or None, optional
            Cell metadata override. If ``None``, uses ``adata.obs``.
    outfile : str, optional
            Path for the downsampled ``.h5ad`` output.
            Default is ``'Group.downsample_1500.h5ad'``.
    downsample : int, optional
            Maximum cells per group. Default is ``1500``.
    """
    adata_path = os.path.expanduser(adata_path)
    outfile = os.path.expanduser(outfile)
    adata = load_adata(adata_path, backed="r")
    if obs_path is not None:
        if isinstance(obs_path, str):
            obs = pd.read_csv(os.path.expanduser(obs_path), sep="\t", index_col=0)
        else:
            obs = obs_path.copy()
        overlapped_cells = list(set(adata.obs_names.tolist()) & set(obs.index.tolist()))
        obs = obs.loc[overlapped_cells]
    else:
        obs = adata.obs.copy()
    keep_cells = (
        obs.loc[obs[groupby].notna()]
        .groupby(groupby)
        .apply(
            lambda x: x.sample(downsample).index.tolist()
            if x.shape[0] > downsample
            else x.index.tolist()
        )
        .sum()
    )
    adata[keep_cells, :].write_h5ad(
        outfile, compression="gzip", convert_strings_to_categoricals=False
    )
    adata.file.close()


def prepare_color_palette(color_dict=None, outpath="palette.xlsx"):
    """
    Generating a .xlsx file including all color palette.

    Parameters
    ----------
    colors : dict
            A dict of dict, keys are categorical terms, values are HEX color code

    Returns
    -------

    """
    outpath = os.path.expanduser(outpath)
    writer = pd.ExcelWriter(outpath)
    for key in color_dict:
        data = pd.DataFrame.from_dict(color_dict[key], orient="index", columns=["Hex"])
        # data.style.background_gradient(cmap='gray_r')
        # data.style.applymap(lambda x:'color:'+x if x.startswith('#') else 'color: white')
        data.to_excel(writer, sheet_name=key, index=True)
        workbook = writer.book
        worksheet = writer.sheets[key]
        colors = data.Hex.tolist()
        for i in range(data.shape[0]):
            color = colors[i]
            f = workbook.add_format(
                {"bold": True, "font_color": "black", "bg_color": color}
            )
            worksheet.write(i + 1, 1, color, f)
        width = 20
        cell_fmt = workbook.add_format(
            {
                "bold": False,
                "font_color": "black",
                # 'bg_color':'green',
                "align": "center",
                "valign": "vcenter",
            }
        )
        # styled = data.style.applymap(lambda val: 'color: %s' % 'red' if val < 0 else 'black').highlight_max()
        worksheet.set_column(0, 1, width, cell_fmt)
    # worksheet.conditional_format(f'A:{last_col}', {'type': 'no_blanks', 'format': cell_fmt})
    writer.close()


def get_color_palette(adata, groupby="Group"):
    """Extract colour palettes from ``adata.obs`` and save to an Excel file.

    Looks for columns named ``color_hex_{level}`` for each of
    ``['Neighborhood', 'Class', 'Subclass', 'Group']`` and writes a
    ``color_palette.xlsx`` file.

    Parameters
    ----------
    adata : AnnData or str
            Annotated data matrix or path to an ``.h5ad`` file.
    groupby : str, optional
            Unused (reserved). Default is ``'Group'``.
    """
    adata = load_adata(adata)
    # color palette
    color_dict = {}
    for col in ["Neighborhood", "Class", "Subclass", "Group"]:
        D = (
            adata.obs.reset_index()
            .loc[:, [col, f"color_hex_{col.lower()}"]]
            .drop_duplicates()
            .set_index(col)[f"color_hex_{col.lower()}"]
            .to_dict()
        )
        color_dict[col] = D
    prepare_color_palette(color_dict=color_dict, outpath="color_palette.xlsx")


def composition(
    obs,
    groupby,
    stratify_col=None,
    composition_col="Region",
    outname=None,
    parent_col=None,
    sort_cols=None,
    adata=None,
    color_palette=None,
):
    """Compute and export cell-type composition across groups to Excel.

    For each stratum (e.g. donor), calculates the proportion of each
    *composition_col* category within each *groupby* group and writes a
    colour-formatted Excel workbook with conditional formatting.

    Parameters
    ----------
    obs : str or DataFrame
            Cell metadata.
    groupby : str
            Column defining the primary grouping (rows in the output).
    stratify_col : str or None, optional
            Column to stratify by (e.g. donor). Each stratum produces a
            separate sheet. If ``None``, only an ``'All'`` sheet is written.
    composition_col : str, optional
            Column whose value distribution is computed within each group.
            Default is ``'Region'``.
    outname : str or None, optional
            Output Excel file name. Defaults to
            ``'{composition_col}_composition.xlsx'``.
    parent_col : str or None, optional
            Higher-level grouping column for hierarchical row annotation.
    sort_cols : list of str or None, optional
            Columns to sort *obs* by before computing order.
    adata : AnnData or str or None, optional
            Used to look up colours from ``adata.uns`` if *color_palette*
            is ``None``.
    color_palette : str, dict, or None, optional
            Colour palette (path to Excel or dict). If ``None`` and *adata*
            is provided, colours are extracted from ``adata.uns``.
    """
    from xlsxwriter.utility import xl_col_to_name

    # for example: Regional Composition Within Each Cell Type (for each cell type in each donor)
    obs = load_obs(obs)
    obs[groupby] = obs[groupby].astype(str)
    if sort_cols is not None:
        obs.sort_values(sort_cols, inplace=True)
    cell_type_order = obs[groupby].unique().tolist()
    if parent_col is not None:
        group2parent = (
            obs.loc[:, [groupby, parent_col]]
            .drop_duplicates()
            .set_index(groupby)[parent_col]
            .to_dict()
        )
    if outname is None:
        outname = f"{composition_col}_composition.xlsx"
    if not outname.endswith(".xlsx"):
        outname = outname + ".xlsx"
    if adata is not None:
        adata = load_adata(adata)
    groups = (
        [parent_col, groupby, composition_col]
        if parent_col is not None
        else [groupby, composition_col]
    )
    if color_palette is not None:
        color_palette = load_color_palette(
            palette=color_palette, adata=adata, groups=groups
        )
    else:
        color_palette = None
    writer = pd.ExcelWriter(outname)
    if stratify_col is not None:  # such as donor
        stratify_order = obs[stratify_col].unique().tolist()
    else:
        stratify_order = []
    for stratify in ["All"] + stratify_order:
        if stratify == "All":
            df = (
                obs.groupby(groupby)[composition_col]
                .value_counts(normalize=True)
                .unstack()
            )
        else:
            df = (
                obs.loc[obs[stratify_col] == stratify]
                .groupby(groupby)[composition_col]
                .value_counts(normalize=True)
                .unstack()
            )
        rows = [ct for ct in cell_type_order if ct in df.index.tolist()]
        df = df.loc[rows]
        if parent_col is not None:
            df.reset_index(inplace=True)
            df.insert(0, parent_col, df[groupby].map(group2parent))
            df.set_index([parent_col, groupby], inplace=True)
        df = df.applymap(lambda x: 100 * x)  # type: ignore
        df.to_excel(writer, sheet_name=stratify)
        workbook = writer.book
        worksheet = writer.sheets[stratify]

        if parent_col is not None:
            parent_values = df.index.get_level_values(0).tolist()
            if color_palette is not None:
                parent_colors = [color_palette[parent_col][ct] for ct in parent_values]
        group_values = df.index.get_level_values(1).tolist()
        regions = df.columns.tolist()
        if color_palette is not None:
            group_colors = [color_palette[groupby][ct] for ct in group_values]
            composition_colors = [color_palette[composition_col][r] for r in regions]
        for i in range(df.shape[0]):
            group_color = group_colors[i] if color_palette is not None else "black"
            # f2 = workbook.add_format({'bold': True, 'font_color': 'black', 'bg_color': group_color})
            f2 = workbook.add_format(
                {"bold": True, "font_color": group_color, "bg_color": "white"}
            )
            if parent_col is not None:
                parent_color = (
                    parent_colors[i] if color_palette is not None else "white"
                )
                f1 = workbook.add_format(
                    {"bold": True, "font_color": "black", "bg_color": parent_color}
                )
                worksheet.write(i + 1, 0, parent_values[i], f1)
                worksheet.write(i + 1, 1, group_values[i], f2)
            else:
                worksheet.write(i + 1, 0, group_values[i], f2)
        for i in range(df.shape[1]):
            composition_color = (
                composition_colors[i] if color_palette is not None else "black"
            )
            f1 = workbook.add_format(
                {"bold": True, "font_color": composition_color, "bg_color": "white"}
            )
            if parent_col is not None:
                worksheet.write(0, i + 2, regions[i], f1)
            else:
                worksheet.write(0, i + 1, regions[i], f1)
        # width = 20
        # cell_fmt = workbook.add_format(
        #     {'bold': False, 'font_color': 'black',
        #      # 'bg_color':'green',
        #      'align': 'center', 'valign': 'vcenter'})
        # worksheet.set_column(0, 1, width, cell_fmt)

        end = str(xl_col_to_name(df.shape[1])) + str(df.shape[0] + 1)
        start = "B2" if parent_col is None else "C2"
        mid_value = 1 / df.shape[1]
        worksheet.conditional_format(
            f"{start}:{end}",
            {
                "type": "3_color_scale",
                "min_color": "#3abf99",
                "min_value": 0,
                "mid_color": "white",
                "mid_value": mid_value,
                "max_color": "#c72228",
                "max_value": 1,
            },
        )
    writer.close()


def taxonomy(
    obs,
    levels=["Neighborhood", "Class", "Subclass", "Group"],
    groupby="Region",
    outfile=None,
):
    """
    _summary_

    Parameters
    ----------
    obs : path or dataframe
            path to obs file (tsv, csv or h5ad) or obs dataframe with
            _description_
    """
    obs = load_obs(obs)
    level = levels[-1]
    D = (
        obs.groupby(level)
        .apply(lambda x: str(x[groupby].value_counts().to_dict()))
        .to_dict()
    )
    obs1 = obs.reset_index().loc[:, levels].drop_duplicates()
    obs1[f"{groupby}.Distrubution"] = obs1[level].map(D)
    obs1.sort_values(levels, inplace=True)  # type: ignore
    if outfile is None:
        outfile = f"{level}.{groupby}_taxonomy.xlsx"
    obs1.to_excel(os.path.expanduser(outfile), index=False)


def get_markers_worker(adata1, obs, level, df_gene, key, outdir, topn, downsample):
    """Run differential expression for one group and return top markers.

    Uses Wilcoxon rank-sum test via scanpy to identify marker genes for
    each cluster at the given *level*. Results are cached to a TSV file.

    Parameters
    ----------
    adata1 : AnnData
            Subset of the annotated data matrix (in memory).
    obs : DataFrame or None
            Cell metadata to assign to ``adata1.obs``. If ``None``, uses
            existing obs.
    level : str
            Column in *obs* defining clusters for differential expression.
    df_gene : DataFrame
            Gene annotation table with a ``'gene_name'`` column used to
            filter marker genes.
    key : str
            Key name for ``sc.tl.rank_genes_groups`` and the output TSV file.
    outdir : str
            Directory to save/load cached marker results.
    topn : int
            Number of top marker genes to return per cluster.
    downsample : int
            Maximum cells per cluster before running the test.

    Returns
    -------
    dict
            ``{cluster_name: [gene1, gene2, ...]}`` with up to *topn* markers.
    """
    import scanpy as sc

    if not os.path.exists(os.path.join(outdir, f"{key}.tsv")):
        use_cells = (
            adata1.obs.groupby(level)
            .apply(
                lambda x: x.sample(downsample).index.tolist()
                if x.shape[0] > downsample
                else x.index.tolist()
            )
            .sum()
        )
        use_adata = adata1[use_cells, :].to_memory()
        if obs is not None:
            use_adata.obs = obs.loc[use_cells, :].copy()
        vc = use_adata.obs[level].value_counts()
        groups = vc[vc >= 3].index.tolist()
        sc.tl.rank_genes_groups(
            use_adata,
            groupby=level,
            groups=groups,
            method="wilcoxon",
            use_raw=False,
            key_added=key,
        )
        # The Wilcoxon rank-sum test is a non-parametric test that compares the ranks of gene expression values between groups. It works best with raw counts because normalization can alter the rank distribution
        markers = sc.get.rank_genes_groups_df(use_adata, group=groups, key=key)
        markers = markers.loc[~markers.names.isna()]
        markers = markers.loc[
            (~markers.logfoldchanges.isna())
            & (markers.scores > 0)
            & (markers.pvals < 0.01)
            & (markers.logfoldchanges > 1)
        ]
        if markers.shape[0] == 0:
            return {}
        markers.sort_values("logfoldchanges", ascending=False, inplace=True)
        markers.to_csv(os.path.join(outdir, f"{key}.tsv"), sep="\t", index=False)
    else:
        markers = pd.read_csv(os.path.join(outdir, f"{key}.tsv"), sep="\t")
    markers.names = markers.names.apply(lambda x: x.split(".")[0])
    markers = (
        markers.loc[markers.names.isin(df_gene["gene_name"].tolist())]
        .groupby("group")
        .apply(lambda x: x.head(topn).names.tolist())
        .to_dict()
    )
    return markers


def get_markers(
    adata_path, levels, gtf, obs=None, outdir="markers", topn=20, downsample=2000
):
    """Identify marker genes at multiple hierarchical levels.

    Iterates over *levels* (e.g. ``['Class', 'Subclass', 'Group']``),
    running Wilcoxon differential expression within each parent group,
    and saves an ``markers.xlsx`` summary.

    Parameters
    ----------
    adata_path : str
            Path to the ``.h5ad`` file.
    levels : list of str
            Hierarchical annotation levels from coarsest to finest.
    gtf : str
            Path to a GTF file for gene name filtering.
    obs : str, DataFrame, or None, optional
            Cell metadata. If ``None``, uses ``adata.obs``.
    outdir : str, optional
            Directory for output files. Default is ``'markers'``.
    topn : int, optional
            Number of top markers per cluster per level. Default is ``20``.
    downsample : int, optional
            Maximum cells per cluster for the test. Default is ``2000``.
    """
    outdir = os.path.expanduser(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    adata = load_adata(adata_path)
    if obs is not None:
        obs = load_obs(obs)
        overlapped_cells = list(set(adata.obs_names.tolist()) & set(obs.index.tolist()))
        obs = obs.loc[overlapped_cells]
        for level in levels:
            adata.obs[level] = adata.obs.index.to_series().map(obs[level].to_dict())
    else:
        obs = adata.obs.copy()
    df_gene = parse_gtf(
        gtf=gtf
    )  # ['chrom','beg','end','gene_name','gene_id','strand','gene_type']
    # DEG
    R = []
    for i in range(len(levels)):
        level = levels[i]
        logger.info(level)
        if i == 0:
            parent = ""
            markers = get_markers_worker(
                adata, obs, level, df_gene, level, outdir, topn, downsample
            )
            for k in markers:  # k is cell type in different levels
                R.append([level, parent, k, markers[k]])  # append topn marker genes
        else:  # i >= 1
            for clusters, df1 in adata.obs.groupby(levels[:i], observed=True):
                if df1.shape[0] == 0:
                    continue
                if df1[level].nunique() < 2:
                    continue
                logger.info(clusters)
                adata1 = adata[df1.index.tolist(), :].to_memory()
                parent = list(clusters)[-1]  # '|'.join(list(clusters))
                key = parent + "|" + level
                # print(key)
                markers = get_markers_worker(
                    adata1, obs, level, df_gene, key, outdir, topn, downsample
                )
                if len(markers) == 0:
                    continue
                for k in markers:
                    R.append([level, parent, k, markers[k]])
    if adata.isbacked:
        adata.file.close()
    data = pd.DataFrame(R, columns=["Level", "Parent", "CellType", "Markers"])
    data.to_excel(
        os.path.join(outdir, "markers.xlsx"), index=False, sheet_name="markers"
    )
