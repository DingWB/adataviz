# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
# import sphinx_rtd_theme,sphinx_sizzle_theme
import os
import sys
import recommonmark # used to parse markdown files
from recommonmark.transform import AutoStructify

# Add project root to sys.path so autodoc can import the package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

project = 'adataviz'
copyright = '2026, Wubin Ding, Amit Klein'
author = 'Wubin Ding, Amit Klein'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "recommonmark",
    "sphinx.ext.napoleon",
    "sphinxcontrib.jquery",
    "sphinx_search.extension",
    "nbsphinx",
]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# select a theme: https://sphinx-themes.org/#themes
# https://sphinx-themes.org/sample-sites/furo/
html_theme = 'furo' # pip install furo
html_static_path = ['_static']


# html_theme = "sphinx_rtd_theme" #sphinx_sizzle_theme
# html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]


source_parsers = {
   '.md': 'recommonmark.parser.CommonMarkParser',
}
source_suffix = ['.rst', '.md']


def setup(app):
    app.add_config_value('recommonmark_config', {
            'auto_toc_tree_section': 'Contents',
            }, True)
    app.add_transform(AutoStructify)