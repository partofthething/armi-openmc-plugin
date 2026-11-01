"""
End-to-end runs of the OpenMC solver at low statistics.

These need the nuclear data and an openmc executable or docker, and take a minute or two.
"""

import numpy as np
import pytest

from armi_openmc import mgxs, runner, verify
from armi_openmc.core import CoreModelBuilder
from armi_openmc.latticePhysics import generateLibrary
from armi_openmc.results import CoreResults
from armi_openmc.tests.conftest import needsSolver

pytestmark = [pytest.mark.openmc, needsSolver]


def test_c5g7SliceMapsPowerAndPins(c5g7, options):
    options.axialSlice = 30.0
    options.pinPowers = True
    model = CoreModelBuilder(c5g7.core, options).build()
    model.model.export_to_model_xml(f"{options.runDir}/model.xml")
    results = CoreResults(model, runner.run(options, options.runDir))
    results.apply(c5g7, 1.0e9)
    assert 1.05 < results.keff < 1.3
    blocks = list(model.blockCells)
    assert sum(b.p.power for b in blocks) == pytest.approx(0.25e9)
    fuel = [b for b in blocks if b.p.openmcPinPowers is not None]
    assert len(fuel) == 4
    for b in fuel:
        assert len(b.p.openmcPinPowers) == 265  # 264 fuel pins and the fission chamber
        assert sum(b.p.openmcPinPowers) == pytest.approx(b.p.power, rel=1e-6)
        assert len(b.p.mgFlux) == len(model.energyBounds) - 1


def test_isotxsReproducesKInfinity(c5g7, options):
    """A homogenized assembly's ISOTXS must give back the k-infinity OpenMC computed for it."""
    options.blockModel = "homogenized"
    options.particles = 5000
    options.batches = 60
    options.inactive = 20
    block = next(b for b in c5g7.core.iterBlocks() if b.getType() == "uo2")
    captured = {}
    original = mgxs.reduce

    def keep(xsModel, statepoint):
        constants, info = original(xsModel, statepoint)
        captured.update(info)
        return constants, info

    mgxs.reduce = keep
    try:
        lib = generateLibrary(c5g7, None, [block], options=options, xsSettings={})
    finally:
        mgxs.reduce = original
    macros = verify.macroscopic(lib, block)
    k = verify.infiniteMediumK(macros, verify.fluxWeightedChi(lib, block, captured["flux"]))
    reference = captured["keff"]
    assert abs(k - reference.nominal_value) < 4 * reference.std_dev + 0.001
    assert np.all(np.isfinite(macros.total))
