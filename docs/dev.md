
# website
```shell
pip install sphinx sphinx-autobuild sphinx-rtd-theme pandoc nbsphinx sphinx_pdj_theme sphinx_sizzle_theme recommonmark readthedocs-sphinx-search

mkdir -p docs && cd docs
sphinx-quickstart
# Separate source and build directories (y/n) [n]: y
# Project name: adataviz

# vim source/conf.py
# add *.rst

#add modules under toctree
sphinx-apidoc -e -o source -f --ext-autodoc --ext-viewcode --ext-githubpages --ext-doctest ../../adataviz #pwd is doc, output (modules.rst) is current dir, source dir is parent dir
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
#还可以使用sphinx-autobuild来自动重载文档。运行下面指令来实现
sphinx-autobuild . _build/html
make latex
```