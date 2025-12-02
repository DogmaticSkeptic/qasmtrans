import os
import sys
from datetime import datetime

project = "QASMTrans"
copyright = f"{datetime.now().year}"
author = "QASMTrans contributors"

extensions = [
    "myst_parser",
    "breathe",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

myst_enable_extensions = ["deflist", "colon_fence"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_static_path = ["_static"]

# Doxygen + Breathe
breathe_projects = {
    "qasmtrans": os.path.abspath(os.path.join("..", "doxygen", "xml")),
}
breathe_default_project = "qasmtrans"

# Ensure source paths are discoverable (if needed later)
sys.path.insert(0, os.path.abspath("../.."))
