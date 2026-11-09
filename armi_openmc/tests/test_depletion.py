"""Depletion with openmc.deplete: which materials burn, names both ways, and the write-back."""

import copy
import os
from contextlib import contextmanager
from types import SimpleNamespace

import pytest
from armi.reactor.flags import Flags

from armi_openmc import depletion, runner
from armi_openmc.core import CoreModelBuilder
from armi_openmc.materials import NuclideMapper
from armi_openmc.options import OpenMCOptions
from armi_openmc.tests.conftest import crossSections, needsData, needsSolver


@contextmanager
def depletableFuel(core):
    """Flag every fuel component DEPLETABLE for the duration (the reactor is shared)."""
    fuel = [c for b in core.iterBlocks() for c in b if c.hasFlags(Flags.FUEL, exact=False)]
    for c in fuel:
        c.p.flags |= Flags.DEPLETABLE
    try:
        yield fuel
    finally:
        for c in fuel:
            c.p.flags &= ~Flags.DEPLETABLE


@needsData
def test_depletableMaterialsAreTheFlaggedComponents(c5g7, options):
    with depletableFuel(c5g7.core) as fuel:
        builder = CoreModelBuilder(c5g7.core, options)
        model = builder.build()
        burnable = depletion.depletableMaterials(model, builder.materials)
    assert {id(c) for _, c, _ in burnable} == {id(c) for c in fuel}
    for mat, c, volume in burnable:
        assert builder.materials.owners[mat.id] is c
        # the quarter is modeled as it is, so the model holds exactly ARMI's share
        assert volume == pytest.approx(c.getVolume())


@needsData
def test_nothingToDepleteIsAnError(c5g7, options):
    options.depletionChain = __file__  # any existing file; the run stops before it is read
    with pytest.raises(ValueError, match="DEPLETABLE"):
        depletion.deplete(c5g7.core, options, 1.0, 1e6, "none")


@needsData
def test_namesBothWays():
    mapper = NuclideMapper(crossSections())
    assert depletion.toOpenMCNames("AM242M", mapper) == [("Am242_m1", 1.0)]
    # isotopes keep their names whether or not the library has data for them
    assert depletion.toOpenMCNames("XE131M", mapper) == [("Xe131_m1", 1.0)]
    oxygen = dict(depletion.toOpenMCNames("O", mapper))
    assert sum(oxygen.values()) == pytest.approx(1.0)
    assert "O16" in oxygen
    assert depletion.toArmiName("Am242_m1") == "AM242M"
    assert depletion.toArmiName("U235") == "U235"


@needsData
def test_writeBackSetsDensitiesAndBurnup(c5g7):
    b = copy.deepcopy(
        next(b for b in c5g7.core.iterBlocks() if any(c.hasFlags(Flags.FUEL) for c in b))
    )
    fuel = next(c for c in b if c.hasFlags(Flags.FUEL))
    b.p.molesHmBOL = b.getHMMoles()
    burned = {nuc: 0.9 * dens for nuc, dens in fuel.getNumberDensities().items()}
    burned = {
        depletion.toOpenMCNames(n, NuclideMapper(crossSections()))[0][0]: d
        for n, d in burned.items()
    }
    burned["Xe135"] = 1e-9
    mat = SimpleNamespace(id=7)
    depletion.writeBack([(mat, fuel, 1.0)], {"7": burned})
    assert fuel.getNumberDensity("XE135") == pytest.approx(1e-9)
    assert b.p.percentBu > 0.0


def test_stepsEndAtEachNode():
    """Node n burns the step that ended there; node 0 of a later cycle decays the outage."""
    o = SimpleNamespace(
        stepLengths=[[10.0, 20.0], [30.0]],
        powerFractions=[[1.0, 0.5], [1.0]],
        cycleLengths=[100.0, 50.0],
        availabilityFactors=[0.8, 1.0],
    )
    fake = SimpleNamespace(o=o, cs={"power": 2e6}, r=None)
    step = depletion.OpenMCDepletionInterface.stepToNode
    assert step(fake, 0, 0) == (0.0, 0.0)
    assert step(fake, 0, 1) == (10.0, 2e6)
    assert step(fake, 0, 2) == (20.0, 1e6)
    assert step(fake, 1, 0) == pytest.approx((20.0, 0.0))
    assert step(fake, 1, 1) == (30.0, 2e6)


@needsData
def test_aProgramRunsInTheContainer(tmp_path):
    options = OpenMCOptions(runner="docker", crossSections=crossSections(), threads=2)
    chainDir = str(tmp_path / "chains")
    cmd, _ = runner.command(
        options,
        tmp_path,
        args=["depletion.json"],
        program=["python", "deplete_step.py"],
        mounts=[chainDir],
        env={"OMP_NUM_THREADS": "2"},
    )
    assert cmd[-3:] == ["python", "deplete_step.py", "depletion.json"]
    assert f"{chainDir}:{chainDir}:ro" in cmd
    assert "OMP_NUM_THREADS=2" in cmd
    native = OpenMCOptions(runner="native", crossSections=crossSections(), python="/opt/py")
    cmd, _ = runner.command(native, tmp_path, program=["python", "deplete_step.py"])
    assert cmd == ["/opt/py", "deplete_step.py"]


def test_autoDepletesWhereOpenMCPythonIs():
    options = OpenMCOptions(runner="auto", python="/no/such/python")
    resolved = depletion.withDepletionRunner(options)
    assert resolved.runner in ("docker", "native")
    if resolved.runner == "docker":
        assert options.runner == "auto"  # the case's options are left alone
    assert depletion.withDepletionRunner(OpenMCOptions(runner="native")).runner == "native"


needsChain = pytest.mark.skipif(
    not os.path.exists(os.environ.get("OPENMC_CHAIN_FILE", "")),
    reason="needs a depletion chain in OPENMC_CHAIN_FILE",
)


def fissile(c):
    return sum(c.getNumberDensity(nuc) for nuc in ("U235", "PU239", "PU241"))


@pytest.mark.openmc
@needsSolver
@needsChain
def test_c5g7SliceBurnsADay(c5g7, options):
    options.axialSlice = 30.0
    options.particles = 500
    options.batches = 15
    options.inactive = 5
    options.depletionIntegrator = "predictor"
    with depletableFuel(c5g7.core) as fuel:
        before = {id(c): fissile(c) for c in fuel}
        saved = {id(c): c.getNumberDensities() for c in fuel}
        try:
            keff = depletion.deplete(c5g7.core, options, 1.0, 1e5, "test")
            assert 0.5 < keff[0] < 1.5
            sliced = [c for c in fuel if c.parent.p.zbottom <= 30.0 < c.parent.p.ztop]
            assert sliced
            for c in sliced:  # only the blocks at the slice are in the model
                if before[id(c)] > 0.0:
                    # not U-235 alone: in MOX the alpha decay of Pu-239 makes more of it
                    assert fissile(c) < before[id(c)]
                    assert c.getNumberDensity("XE135") > 0.0
        finally:
            for c in fuel:
                c.setNumberDensities(saved[id(c)])
