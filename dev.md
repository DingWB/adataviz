
# Update website
```shell
pip install sphinx sphinx-autobuild sphinx-rtd-theme pandoc nbsphinx sphinx_pdj_theme sphinx_sizzle_theme recommonmark readthedocs-sphinx-search
conda install conda-forge::pandoc

mkdir -p docs && cd docs
sphinx-quickstart
# Separate source and build directories (y/n) [n]: y
# Project name: adataviz

# vim source/conf.py
# add *.rst

# cd docs
# vim index.html: <meta http-equiv="refresh" content="0; url=./build/html/index.html" />
cd docs
rm -rf build
ln -s ~/Projects/Github/adataviz/notebooks source/notebooks
sphinx-apidoc -e -o source -f ../../adataviz
make html
rm -rf source/notebooks
cd ..
ls
ls docs

vim .nojekyll #create empty file
```

# Hooks
```shell
# install pixi
curl -fsSL https://pixi.sh/install.sh | bash
export PATH="$HOME/.pixi/bin:$PATH" && which pixi
pixi run hooks-install
```

# Upload onto conda
## Step1: publish onto pypi
```shell
curl -sL https://pypi.io/packages/source/a/adataviz/adataviz-0.5.tar.gz | sha256sum
# put sha256 in meta.yaml for the first time
```

## Step2
```shell
# 1. Fork https://github.com/conda-forge/staged-recipes
git clone https://github.com/DingWB/staged-recipes.git
cd staged-recipes
git checkout -b adataviz # create a new branch

# 2. add recipe
recipes/adataviz/meta.yaml

git add recipes/adataviz/meta.yaml
git commit -m "Add adataviz recipe"
git push origin adataviz
# PR to bioconda/bioconda-recipes:master, and finally:
conda install -c conda-forge adataviz
# @conda-forge/help-python: AFAICT, this is ready for review.
```

Run grayskull to generate new version (recipe.yaml)
```shell
pip install grayskull
cd ~/Projects/Github/staged-recipes/recipes
grayskull pypi --strict-conda-forge adataviz

cd /home/x-wding2/Projects/Github/staged-recipes
pixi run lint #recipes/adataviz/recipe.yaml
# or conda-smithy
conda smithy recipe-lint recipes/adataviz
```