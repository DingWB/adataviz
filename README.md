# adataviz
Functions and tools to visualize adata

## **Installation**
----------------------
1. **Install using pip**:
```shell
pip install adataviz

#upgrade from older version
pip install --upgrade adataviz
```

2. **Install the developmental version directly from github**:
```shell
pip install git+https://github.com/DingWB/adataviz
# reinstall
pip uninstall -y adataviz && pip install git+https://github.com/DingWB/adataviz

```
OR
```shell
git clone https://github.com/DingWB/adataviz
cd adataviz
python setup.py install
```

## Command Line Tools
### to_pseudobulk
```shell
# merge single cells raw counts to psuedobulk (sum up raw counts) and run normalization (logCPM or logTPM)
adataviz tool scrna2pseudobulk  HMBA.Group.downsample_1500.h5ad --groupby="Subclass" --downsample=2000 --use_raw=True -m RNA --n_jobs 16 --normalization CPM -s ~/Projects/BICAN/adata/HMBA_v2/Pseudobulk.Subclass.h5ad

# stat pseudobulk: calculate min,q25, q50, q75, max, mean and std
adataviz tool stat_pseudobulk  HMBA.Group.downsample_1500.h5ad --groupby="Subclass" --downsample=2000 --use_raw=False -m RNA --n_jobs 16 -s ~/Projects/BICAN/adata/HMBA_v2/Pseudobulk.Subclass.h5ad
```