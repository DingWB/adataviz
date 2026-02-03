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
html_css_files = [
    'custom.css',
]
# html_js_files = [
#     'plotly_scroll.js',
# ]

html_theme_options = {
    "navigation_with_keys": True,
    "top_of_page_buttons": ["view", "edit"],
    "announcement": "<em>Documentation website</em> is online!",
    # "analytics_id": "G-VRB2NBWG05",
    "collapse_navigation": False,
    "globaltoc_collapse": False,
    "globaltoc_maxdepth": 3,
    "collapse_navigation": False,
    "display_version": True,
    # "sidebarwidth": 200,  # sidebarwidth
    "navigation_depth": 6,
}

# Use MathJax v2 site-wide to match notebook fragments (avoid v2/v3 conflict)
# nbsphinx/nbconvert output embeds MathJax v2; forcing Sphinx to use the
# same path avoids loading multiple MathJax versions which can cause
# "Cannot read properties of undefined (reading 'Startup')" errors.
mathjax_path = (
    "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_SVG"
)


# html_theme = "sphinx_rtd_theme" #sphinx_sizzle_theme
# html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# html_sidebars = {
#     "**": [
#         "relations.html",  # needs 'show_related': True theme option to display
#         "searchbox.html",
#     ]
# }

html_context = {
    "display_github": True,
    "github_user": "DingWB",
    "github_repo": "adataviz",
    "github_version": "main/docs/source/",
}

source_parsers = {
   '.md': 'recommonmark.parser.CommonMarkParser',
}
source_suffix = ['.rst', '.md']


def setup(app):
    app.add_config_value('recommonmark_config', {
            'auto_toc_tree_section': 'Contents',
            }, True)
    app.add_transform(AutoStructify)




