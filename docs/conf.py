import os
import sys

project = "ClassyMC"
copyright = "ClassyMC Developers"
author = "ClassyMC Developers"
release = "1.0"

extensions = [
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "navigation_depth": 4,
    "titles_only": False,
    "collapse_navigation": False,
    "sticky_navigation": True,
}

html_css_files = ["custom.css"]

autosectionlabel_prefix_document = True
