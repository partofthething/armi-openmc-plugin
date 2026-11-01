"""``python -m armi_openmc <command> ...``: the ARMI command line with the OpenMC plugin loaded."""

import sys

import armi


def main():
    from armi_openmc.app import OpenMCApp

    if not armi.isConfigured():
        armi.configure(OpenMCApp())
    from armi.cli import ArmiCLI

    sys.exit(ArmiCLI().run())


if __name__ == "__main__":
    main()
