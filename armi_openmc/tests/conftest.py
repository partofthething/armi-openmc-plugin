"""Shared fixtures: the example C5G7 and FFTF reactors, and the nuclear data."""

import os
import pathlib
import shutil

import armi
import pytest

if not armi.isConfigured():
    from armi_openmc.app import OpenMCApp

    armi.configure(OpenMCApp())

ROOT = pathlib.Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"
DEFAULT_DATA = "/mnt/storage/nuclear_data/openmc/endfb-viii.1-hdf5/cross_sections.xml"


def crossSections():
    path = os.environ.get("OPENMC_CROSS_SECTIONS", DEFAULT_DATA)
    return path if os.path.exists(path) else None


needsData = pytest.mark.skipif(crossSections() is None, reason="no OpenMC nuclear data")
needsSolver = pytest.mark.skipif(
    crossSections() is None or not (shutil.which("openmc") or shutil.which("docker")),
    reason="needs nuclear data and an openmc executable or docker",
)


def loadReactor(directory, settingsFile):
    from armi import init

    previous = os.getcwd()
    os.chdir(EXAMPLES / directory)
    try:
        return init(fName=settingsFile).r
    finally:
        os.chdir(previous)


@pytest.fixture(scope="session")
def c5g7():
    return loadReactor("c5g7", "c5g7-settings.yaml")


@pytest.fixture(scope="session")
def fftf():
    return loadReactor("fftf", "FFTF.yaml")


@pytest.fixture
def options(tmp_path):
    from armi_openmc.options import OpenMCOptions

    return OpenMCOptions(
        crossSections=crossSections() or "",
        particles=2000,
        batches=30,
        inactive=10,
        runDir=str(tmp_path),
    )
