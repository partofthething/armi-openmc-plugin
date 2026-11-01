"""
OpenMC cross sections for the mixture behind ARMI's MC2-3 ISOTXS fixture, compared with it.

ARMI's test fixture ``armi/nuclearDataIO/tests/fixtures/mc2v3-AA.isotxs`` was made by MC2-3
(ENDF/B-VII.0, ANL33) from the infinite homogeneous fast mixture in
``armi/nuclearDataIO/tests/library-file-generation/mc2v3-AA.inp``. This builds the same mixture as
an ARMI block, makes its ISOTXS with the OpenMC plugin and compares the two.

Differences to expect, apart from statistics: ENDF/B-VII.0 vs VIII.1, MC2-3's critical
buckling search (it adds leakage, hardening the spectrum; OpenMC runs the infinite medium), and
one material temperature here (873 K) where MC2-3 had 743 K iron and 738 K sodium.

    python compare_mc2v3.py [particles] [batches]
"""

import os
import sys

import armi

if not armi.isConfigured():
    from armi_openmc.app import OpenMCApp

    armi.configure(OpenMCApp())

from armi.reactor.blocks import HexBlock  # noqa: E402
from armi.reactor.components import Hexagon  # noqa: E402

from armi_openmc import cli, mgxs  # noqa: E402
from armi_openmc.latticePhysics import generateLibrary  # noqa: E402
from armi_openmc.options import OpenMCOptions  # noqa: E402

FIXTURE = os.path.join(
    os.path.dirname(armi.__file__), "nuclearDataIO", "tests", "fixtures", "mc2v3-AA.isotxs"
)
DENSITIES = {
    "U235": 1.0e-3,
    "U238": 1.0e-3,
    "PU239": 1.0e-4,
    "FE54": 1.0e-3,
    "FE56": 1.0e-2,
    "FE57": 1.0e-4,
    "FE58": 1.0e-5,
    "NA23": 1.0e-3,
    "ZR90": 1.0e-3,
    "ZR91": 1.0e-4,
    "ZR92": 1.0e-4,
    "ZR94": 1.0e-4,
    "ZR96": 1.0e-5,
}
TEMPERATURE_C = 873.0 - 273.15


def mixtureBlock():
    """A one-component hexagonal block of the MC2-3 mixture, cross section ID AA."""
    c = Hexagon(
        "mixture", "Void", Tinput=TEMPERATURE_C, Thot=TEMPERATURE_C, ip=0.0, op=16.0, mult=1
    )
    b = HexBlock("mixture")
    b.add(c)
    b.setHeight(10.0)
    c.setNumberDensities(DENSITIES)
    b.setType("mixture")
    b.p.xsType = "A"
    b.p.envGroup = "A"
    return b


def main(particles=20000, batches=150):
    options = OpenMCOptions(
        crossSections=os.environ.get(
            "OPENMC_CROSS_SECTIONS",
            "/mnt/storage/nuclear_data/openmc/endfb-viii.1-hdf5/cross_sections.xml",
        ),
        particles=int(particles),
        batches=int(batches),
        inactive=30,
        blockModel="homogenized",
        groupStructure="ANL33",
        runDir="openmc-runs",
    )
    block = mixtureBlock()
    lib = generateLibrary(None, None, [block], options=options, xsSettings={})
    mgxs.writeLibrary(lib, "ISOAA-openmc")
    print("OpenMC (ENDF/B-VIII.1) against MC2-3 (ENDF/B-VII.0), collapsed with the OpenMC flux:")
    cli.compareLibraries(FIXTURE, "ISOAA-openmc")


if __name__ == "__main__":
    main(*sys.argv[1:])
