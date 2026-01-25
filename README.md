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
adataviz tool to_pseudobulk  HMBA.Group.downsample_1500.h5ad --groupby="Subclass" --downsample=2000 --use_raw=True -m RNA --n_jobs 16 --normalization CPM -s ~/Projects/BICAN/adata/HMBA_v2/Pseudobulk.Subclass.h5ad
```