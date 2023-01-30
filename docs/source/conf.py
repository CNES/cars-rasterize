#!/usr/bin/env python
# coding: utf8
#
# Configuration file for the Sphinx documentation builder.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# pylint: skip-file
# flake8: noqa
# type: ignore

import os
import sys

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath("../.."))

# Extend Recursion limit for RecursionError in big files (bug astroid)
sys.setrecursionlimit(8 * sys.getrecursionlimit())

# -- Project information -----------------------------------------------------

# General information about the project.
project = "cars-rasterize"
copyright = "2023, Centre National d'Etudes Spatiales (CNES)"
author = "CNES"

# The full version, including alpha/beta/rc tags
from pkg_resources import get_distribution

try:
    version = get_distribution("cars_rasterize").version
    release = version
except Exception as error:
    print("WARNING: cannot find cars_rasterize version")
    version = "Unknown"
    release = version

# The master toctree document.
master_doc = "index"

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.intersphinx",  # other projects automatic links to doc
    "sphinx.ext.imgmath",  # Add rst math capabilities with :math:
    "sphinx.ext.autodoc",  # apidoc automatic generation
    "sphinx.ext.viewcode",  # viewcode in automatic apidoc
    "autoapi.extension",
    "sphinx_tabs.tabs",
]

# imgmath configuration
imgmath_embed = True

# Autoapi configuration
autoapi_dirs = ["../../cars_rasterize"]
autoapi_root = "api_reference"
autoapi_keep_files = True
autoapi_options = [
    "members",
    "undoc-members",
    "private-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
]

# add autohint
autodoc_typehints = "description"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Title
html_title = "cars-rasterize Documentation"
html_short_title = "cars-rasterize Documentation"

# Logo
# html_logo = "images/picto_transparent_mini.png"

# Favicon
# html_favicon = "images/favicon_noname.ico"

# Theme options
html_theme_options = {
    "logo_only": True,
    "navigation_depth": 3,
}

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = ["css/my_custom.css"]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    "papersize": "letterpaper",
    # The font size ('10pt', '11pt' or '12pt').
    "pointsize": "10pt",
    # Additional stuff for the LaTeX preamble.
    "preamble": "",
    # Latex figure (float) alignment
    "figure_align": "htbp",
}
numfig = True

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        master_doc,
        "cars_rasterize.tex",
        "cars-rasterize documentation",
        "TODO",
        "manual",
    ),
]
