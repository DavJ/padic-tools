import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
project = "padic-tools"
author = "David Jaroš"
copyright = "2026, David Jaroš"
release = "0.1.0"
extensions = [
    "sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.doctest",
    "sphinx.ext.viewcode", "sphinxcontrib.katex",
]
html_theme = "alabaster"
html_theme_options = {"description": "p-adic analysis · exact algebra · explicit bridges", "fixed_sidebar": True}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
exclude_patterns = ["_build"]
autodoc_member_order = "bysource"
autodoc_typehints = "description"
# Local JavaScript plus native MathML: equations need no CDN or webfonts.
katex_css_path = "katex-native.css"
katex_options = 'output: "mathml", throwOnError: true'
katex_prerender = False
doctest_global_setup = "from fractions import Fraction; from padic_tools import *"
