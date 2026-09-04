from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

project = "SpPrAn"
author = "Sergio Zamora-Erazo"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_rtd_theme",
]

autodoc_mock_imports = [
    "anndata", "matplotlib", "networkx", "numpy", "pandas", "pyarrow",
    "scanpy", "scipy", "seaborn", "squidpy", "statsmodels", "yaml",
]

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

napoleon_google_docstring = False
napoleon_numpy_docstring = True
templates_path = ["_templates"]
exclude_patterns = []
html_theme = "sphinx_rtd_theme"
intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}
