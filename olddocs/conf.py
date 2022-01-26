#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Largely adapted from Numpy's conf.py

import sys, os, re

sys.path.insert(0, os.path.abspath("../src"))

# -- General configuration -----------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.graphviz",
    "sphinx.ext.ifconfig",
    "matplotlib.sphinxext.plot_directive",
    "numpydoc",
]

autosummary_generate = True

extensions.append("sphinx.ext.imgmath")
imgmath_image_format = "svg"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# autoclass_content = 'both'


# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "Alhambra"
copyright = "2021, Constantine Evans"
author = "Constantine Evans"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#

import alhambra

version = re.sub(r"(\d+\.\d+)\.\d+(.*)", r"\1\2", alhambra.__version__)
version = re.sub(r"(\.dev\d+).*?$", r"\1", version)
# The full version, including alpha/beta/rc tags.
release = alhambra.__version__


add_function_parentheses = False

default_role = "autolink"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

pygments_style = "sphinx"


# -- Options for HTML output ----------------------------------------------

# themedir = os.path.join(os.path.dirname(__file__), "scipy-sphinx-theme", "_theme")
# if not os.path.isdir(themedir):
#    raise RuntimeError(
#        "Get the scipy-sphinx-theme first, "
#        "via git submodule init && git submodule update"
#    )

# html_theme = "scipy"
# html_theme_path = [themedir]


# Default build
# html_theme_options = {
#     "edit_link": False,
#     "sidebar": "left",
#     "scipy_org_logo": False,
#     "rootlinks": [],
# }
html_sidebars = {}

html_title = "%s v%s Manual" % (project, version)
html_last_updated_fmt = "%b %d, %Y"
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "alhambradoc"


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    "papersize": "a4paper",
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        master_doc,
        "alhambra.tex",
        "Alhambra Documentation",
        "Constantine Evans",
        "manual",
    ),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "alhambra", "Alhambra Documentation", [author], 1)]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "alhambra",
        "alhambra Documentation",
        author,
        "alhambra",
        "One line description of project.",
        "Miscellaneous",
    ),
]

numpydoc_show_inherited_class_members = False
# numpydoc_class_members_toctree = False
