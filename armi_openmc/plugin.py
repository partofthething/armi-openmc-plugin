"""The ARMI plugin: settings, parameters, interfaces and the input editor's OpenMC section."""

from armi import plugins


class OpenMCPlugin(plugins.ArmiPlugin):
    """OpenMC as ARMI's global flux solver and lattice physics code."""

    @staticmethod
    @plugins.HOOKIMPL
    def exposeInterfaces(cs):
        from armi_openmc import globalFlux, latticePhysics

        info = []
        for mod in (latticePhysics, globalFlux):
            info += plugins.collectInterfaceDescriptions(mod, cs)
        return info

    @staticmethod
    @plugins.HOOKIMPL
    def defineSettings():
        from armi_openmc import settings

        return settings.defineSettings()

    @staticmethod
    @plugins.HOOKIMPL
    def defineParameters():
        from armi_openmc import parameters

        return parameters.getParameterDefinitions()

    @staticmethod
    @plugins.HOOKIMPL
    def defineEntryPoints():
        from armi_openmc import cli

        return [cli.BuildOpenMCModel, cli.GenerateXS, cli.CompareXS]

    @staticmethod
    def defineSettingEditors(cs):
        """Not a hook: the input editor asks plugins for this directly."""
        from armi_openmc import settings

        return settings.defineSettingEditors(cs)
