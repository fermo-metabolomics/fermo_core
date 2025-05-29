from importlib import metadata

project = "fermo_core"
copyright = "2025, Mitja M. Zdouc"
author = "Mitja M. Zdouc"
version = metadata.version("fermo_core")

extensions = ["sphinx_rtd_theme", "autoapi.extension"]

html_theme = "sphinx_rtd_theme"

autoapi_dirs = ["../../fermo_core"]

templates_path = ["_templates"]
html_static_path = ["_static"]

rst_epilog = f"""
.. |version| replace:: {version}
"""
