
# website
```shell
pip install sphinx sphinx-autobuild sphinx-rtd-theme pandoc nbsphinx sphinx_pdj_theme sphinx_sizzle_theme recommonmark readthedocs-sphinx-search
conda install conda-forge::pandoc

mkdir -p docs && cd docs
sphinx-quickstart
# Separate source and build directories (y/n) [n]: y
# Project name: adataviz

# vim source/conf.py
# add *.rst

#add modules under toctree
# cd docs
# sphinx-apidoc -e -o source -f --ext-autodoc --ext-viewcode --ext-githubpages --ext-doctest ../../adataviz #pwd is doc, output (modules.rst) is current dir, source dir is parent dir

# under docs
# vim index.html: <meta http-equiv="refresh" content="0; url=./build/html/index.html" />
cd docs
rm -rf build
ln -s ~/Projects/Github/adataviz/notebooks/ source/notebooks
sphinx-apidoc -e -o source -f ../../adataviz
make html
rm -rf source/notebooks
cd ..
ls

vim .nojekyll #create empty file
make html # make clean && make html
#or run sphinx-autobuild to auto-reload
sphinx-autobuild . _build/html
make latex
```