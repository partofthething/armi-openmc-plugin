"""Sphinx configuration for the armi-openmc docs."""

import pathlib
import sys

import armi

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

if not armi.isConfigured():
    from armi_openmc.app import OpenMCApp

    armi.configure(OpenMCApp())

import armi_openmc  # noqa: E402

project = "armi-openmc"
author = "Nick Touran"
copyright = "2026, Nick Touran"
release = armi_openmc.__version__
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]
autodoc_member_order = "bysource"
autodoc_inherit_docstrings = False  # ARMI's hook docstrings use sphinx-needs directives
autodoc_default_options = {"members": True, "show-inheritance": True}
napoleon_numpy_docstring = True
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "openmc": ("https://docs.openmc.org/en/stable", None),
    "armi": ("https://partofthething.github.io/armi", None),
}
exclude_patterns = ["_build", "make_figures.py"]
nitpicky = False

html_theme = "sphinx_book_theme"
html_title = "ARMI + OpenMC"
html_static_path = ["_static"]
html_last_updated_fmt = "%Y-%m-%d"
html_baseurl = "https://partofthething.github.io/armi-openmc-plugin/"
html_theme_options = {
    "repository_url": "https://github.com/partofthething/armi-openmc-plugin",
    "repository_branch": "main",
    "path_to_docs": "doc",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": True,
    "show_toc_level": 2,
    "home_page_in_toc": True,
    "use_download_button": False,
    "use_fullscreen_button": False,
}


def _settingsTable():
    """The plugin's settings as a list-table, written where ``settings.rst`` includes it."""
    from armi.settings import setting

    from armi_openmc import settings

    lines = [
        ".. list-table::",
        "   :header-rows: 1",
        "   :widths: 22 14 64",
        "",
        "   * - Setting",
        "     - Default",
        "     - What it does",
    ]
    for item in settings.defineSettings():
        if not isinstance(item, setting.Setting):
            continue  # an option the plugin adds to a framework setting
        default = item.default if item.default not in ("", []) else "(blank)"
        options = getattr(item, "options", None)
        text = " ".join(str(item.description).split())
        if options:
            text += " Options: " + ", ".join(f"``{o}``" for o in options) + "."
        lines += [f"   * - ``{item.name}``", f"     - ``{default}``", f"     - {text}"]
    out = pathlib.Path(__file__).parent / "_generated"
    out.mkdir(exist_ok=True)
    (out / "settings-table.rst").write_text("\n".join(lines) + "\n")


_settingsTable()
