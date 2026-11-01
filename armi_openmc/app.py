"""An ARMI application with the OpenMC plugin registered on top of the built-in plugins."""

from armi.apps import App

from armi_openmc.plugin import OpenMCPlugin


class OpenMCApp(App):
    name = "armi-openmc"

    def __init__(self):
        App.__init__(self)
        self._pm.register(OpenMCPlugin)

    @property
    def splashText(self):
        return """
+===============================================+
|   ARMI + OpenMC                               |
|   Monte Carlo for ARMI reactor models         |
+===============================================+
"""
