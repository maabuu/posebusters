"""Configuration file for the Sphinx documentation builder."""

# -- Project information -----------------------------------------------------

project = "PoseBusters"
copyright = "2023, Martin Buttenschoen"
author = "Martin Buttenschoen"

# get version and release from library. ref https://peps.python.org/pep-0440/
version, release = "major.minor.micro", "major.minor"
with open("../../posebusters/__init__.py") as init_file:
    for line in init_file:
        if line.startswith("__version__"):
            version = line.split("=")[1].strip().strip('"')
            release = ".".join(version.split(".", 2)[:2])
            break

# -- General configuration ---------------------------------------------------

extensions = [
    "hoverxref.extension",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    # "sphinx.ext.autosectionlabel",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "nbsphinx",
    "sphinxcontrib.programoutput",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]
exclude_patterns = ["build", "_build", "Thumbs.db", ".DS_Store"]

language = "en"

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "collapse_navigation": False,
    "display_version": True,
    "prev_next_buttons_location": "both",
    "sticky_navigation": True,
    # 'analytics_id': 'UA-45051049-3',
    # "navigation_depth": 2,
    "includehidden": False,
}
html_static_path = ["_static"]
html_logo = "./_static/logo_square.png"
htmlhelp_basename = "posebusters_doc"
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

# -- Options for programoutput extension -------------------------------------
programoutput_prompt_template = ">>> {command}\n{output}"
