"""Configuration file for the Sphinx documentation builder."""

# -- Project information -----------------------------------------------------

project = "MolBuster"
copyright = "2023, Martin Buttenschoen"
author = "Martin Buttenschoen"

release = "0.1"
version = "0.1.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "hoverxref.extension",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "nbsphinx",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]
exclude_patterns = ["build"]

language = "en"

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "collapse_navigation": False,
    "display_version": False,
    "prev_next_buttons_location": "both",
    "sticky_navigation": True,
    # 'analytics_id': 'UA-45051049-3',
    # 'navigation_depth': 2,
}
html_static_path = ["_static"]
html_logo = "./_static/logo_square.png"
htmlhelp_basename = "molbuster_doc"
html_last_updated_fmt = "%d %B %Y"

# -- Options for EPUB output -------------------------------------------------

epub_show_urls = "footnote"

# -- Options for todo extension ----------------------------------------------

todo_include_todos = True

# -- Options for napoleon extension ------------------------------------------

# https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# -- Options for coverage extension ------------------------------------------
coverage_show_missing_items = True

# -- Options for hoverxref extension -----------------------------------------
hoverxref_auto_ref = True
hoverxref_role_types = {
    "class": "tooltip",
    "command": "tooltip",
    "confval": "tooltip",
    "hoverxref": "tooltip",
    "mod": "tooltip",
    "ref": "tooltip",
    "reqmeta": "tooltip",
    "setting": "tooltip",
    "signal": "tooltip",
}
hoverxref_roles = ["command", "reqmeta", "setting", "signal"]
