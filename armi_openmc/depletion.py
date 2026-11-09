"""
Burning the reactor with ``openmc.deplete`` (``openmcDepletion: true``).

ARMI calls its depletion interfaces at every time node, before the flux solve. At node ``n`` of a
cycle this interface depletes the reactor over the step that ended there
(``stepLengths[cycle][n-1]`` full-power days, at that step's power):

1. the whole-core model is built from the reactor as it stands, exactly as the global flux solver
   builds it, and every material made of one depletable component (``isDepletable``) is marked
   depletable, with its volume in the model;
2. ``depletion.json`` gives those materials' full inventory (nuclides without cross section data
   included, which OpenMC carries through decay only), the chain, the integrator, the step and the
   power, and ``deplete_step.py`` runs one ``openmc.deplete`` step in the OpenMC container (or a
   native install) with them;
3. the number densities at the end of the step go back onto the components, every nuclide of the
   chain that ARMI knows, and each fuel block's ``percentBu`` is updated from its heavy metal.

At node 0 of every cycle after the first, the downtime of the cycle before
(``cycleLength * (1 - availabilityFactor)`` days) is decayed at zero power, with no transport.

The interface does not set keff or the power distribution: pair it with ``neutronicsKernel: OpenMC``
(or another flux solver), which runs after it at every node on the depleted reactor.

Only components that are their own OpenMC material can be depleted, i.e. those of heterogeneous
blocks. A depletable component in a homogenized block keeps its number densities, with a warning.
"""

import dataclasses
import functools
import json
import os
import shutil
import subprocess

import openmc.data
from armi import interfaces, runLog
from armi.nucDirectory import nuclideBases
from armi.physics.neutronics.isotopicDepletion.isotopicDepletionInterface import isDepletable

from armi_openmc import runner
from armi_openmc.core import CoreModelBuilder
from armi_openmc.options import OpenMCOptions
from armi_openmc.settings import CONF_DEPLETION, INTEGRATORS

ORDER = interfaces.STACK_ORDER.DEPLETION
#: Nuclides whose fission and capture are tallied directly in the flux reaction rate mode: the
#: major actinides, whose resonance self-shielding a group flux collapse would miss.
DIRECT_NUCLIDES = [
    "U234",
    "U235",
    "U236",
    "U238",
    "Np237",
    "Pu238",
    "Pu239",
    "Pu240",
    "Pu241",
    "Pu242",
    "Am241",
    "Th232",
    "U233",
]
DIRECT_REACTIONS = ["fission", "(n,gamma)"]
SCRIPT = os.path.join(os.path.dirname(__file__), "resources", "deplete_step.py")
SPEC_FILE = "depletion.json"
RESULT_FILE = "depletion-result.json"


def describeInterfaces(cs):
    if cs[CONF_DEPLETION]:
        return (OpenMCDepletionInterface, {})
    return None


class OpenMCDepletionInterface(interfaces.Interface):
    """Depletes the depletable components with openmc.deplete at every time node."""

    name = "openmcDepletion"
    purpose = "depletion"

    def interactBOL(self):
        opts = OpenMCOptions.fromSettings(self.cs)
        runLog.info(
            f"OpenMC depletion: chain {opts.depletionChainPath}, {opts.depletionIntegrator} "
            f"integrator, runner {opts.runner}"
        )

    def interactEveryNode(self, cycle, node):
        days, power = self.stepToNode(cycle, node)
        if days <= 0.0:
            return
        what = f"{days:.2f} days at {power / 1e6:.4g} MW" if power else f"{days:.2f} days of decay"
        runLog.important(f"OpenMC depletion: burning to cycle {cycle} node {node}, {what}")
        keff = deplete(
            self.r.core, OpenMCOptions.fromSettings(self.cs), days, power, f"c{cycle}n{node}"
        )
        if power and keff is not None:
            runLog.info(
                f"OpenMC depletion: keff at the start of the step {keff[0]:.5f} +/- {keff[1]:.5f}"
            )

    def stepToNode(self, cycle, node):
        """Full-power days and power (W) of the step ending at this node (zero for downtime)."""
        o = self.o
        if node > 0:
            basicPower = self.cs["power"] or (
                self.cs["powerDensity"] * self.r.core.getHMMass() * self.r.core.powerMultiplier
            )
            return o.stepLengths[cycle][node - 1], o.powerFractions[cycle][node - 1] * basicPower
        if cycle == 0:
            return 0.0, 0.0
        previous = cycle - 1
        return o.cycleLengths[previous] * (1.0 - o.availabilityFactors[previous]), 0.0


def deplete(core, options: OpenMCOptions, days, power, label):
    """
    Deplete ``core`` for ``days`` at ``power`` W (zero: decay only) and write the result back.

    Returns keff and its standard deviation at the start of the step (None for decay).
    """
    if options.multigroup:
        raise ValueError("openmcDepletion needs continuous-energy OpenMC, not the multigroup mode")
    workdir = os.path.join(options.runDir, f"deplete-{label}")
    os.makedirs(workdir, exist_ok=True)
    for stale in (RESULT_FILE, "depletion_results.h5"):
        if os.path.exists(os.path.join(workdir, stale)):
            os.remove(os.path.join(workdir, stale))
    builder = CoreModelBuilder(core, options)
    coreModel = builder.build()
    burnable = depletableMaterials(coreModel, builder.materials)
    if not burnable:
        raise ValueError(
            "openmcDepletion is on but no component can be depleted: flag the fuel components "
            "DEPLETABLE (`burn: true` in the blueprint's nuclide flags, or `depletable` in their "
            "flags) and model their blocks heterogeneously."
        )
    coreModel.model.export_to_model_xml(os.path.join(workdir, "model.xml"))
    chain = options.depletionChainPath
    spec = {
        "chain": chain,
        "integrator": INTEGRATORS[options.depletionIntegrator],
        "days": float(days),
        "power": float(power) * coreModel.domainFraction,
        "reactionRates": {
            "mode": options.depletionReactionRates,
            "groups": options.depletionGroups,
            "nuclides": DIRECT_NUCLIDES,
            "reactions": DIRECT_REACTIONS,
        },
        "materials": {
            str(mat.id): {"volume": volume, "nuclides": inventory(c, builder.mapper)}
            for mat, c, volume in burnable
        },
    }
    with open(os.path.join(workdir, SPEC_FILE), "w") as f:
        json.dump(spec, f)
    shutil.copy(SCRIPT, os.path.join(workdir, os.path.basename(SCRIPT)))
    threads = options.threads or runner.physicalCores()
    runner.execute(
        withDepletionRunner(options),
        workdir,
        args=[SPEC_FILE],
        logName="deplete.out",
        program=["python", os.path.basename(SCRIPT)],
        mounts=[os.path.dirname(chain)],
        env={"OMP_NUM_THREADS": str(threads)},
    )
    with open(os.path.join(workdir, RESULT_FILE)) as f:
        result = json.load(f)
    writeBack(burnable, result["materials"])
    return tuple(result["keff"]) if power else None


def withDepletionRunner(options: OpenMCOptions) -> OpenMCOptions:
    """
    The options with ``auto`` resolved for depletion.

    Depletion needs OpenMC's Python API and shared library, not just the ``openmc`` executable, so
    ``auto`` runs natively only when ``openmcPython`` can import ``openmc.lib``, and otherwise in
    the container.
    """
    if options.runner != "auto":
        return options
    native = _canDeplete(options.python) or not shutil.which("docker")
    return dataclasses.replace(options, runner="native" if native else "docker")


@functools.lru_cache(maxsize=8)
def _canDeplete(python):
    if not shutil.which(python):
        return False
    check = subprocess.run([python, "-c", "import openmc.lib"], capture_output=True)
    return check.returncode == 0


def depletableMaterials(coreModel, factory):
    """(material, component, model volume) of every material made of one depletable component."""
    owners = factory.owners
    out, skipped = [], set()
    for mat in coreModel.model.materials:
        c = owners.get(mat.id)
        if c is None or not isDepletable(c):
            continue
        b = c.parent
        share = coreModel.copies.get(b, 1) * (1.0 if coreModel.cut else b.getSymmetryFactor())
        volume = c.getVolume() * share
        if volume > 0.0:
            out.append((mat, c, volume))
    depleted = {id(c) for _, c, _ in out}
    for b in coreModel.blockCells:
        for c in b:
            if isDepletable(c) and id(c) not in depleted:
                skipped.add(b.getType())
    if skipped:
        runLog.warning(
            f"OpenMC depletion: depletable components of {sorted(skipped)} blocks are not their "
            "own material (homogenized blocks?) and are not depleted.",
            single=True,
        )
    return out


def inventory(component, mapper):
    """A component's number densities by OpenMC nuclide name, in the library or not."""
    out = {}
    for name, density in component.getNumberDensities().items():
        if density <= 0.0:
            continue
        for nuc, frac in toOpenMCNames(name, mapper):
            out[nuc] = out.get(nuc, 0.0) + density * frac
    return out


def toOpenMCNames(armiName, mapper):
    """OpenMC (GNDS) names and fractions of one ARMI nuclide, in the library or not."""
    base = nuclideBases.byName.get(armiName)
    if base is None or isinstance(
        base, (nuclideBases.LumpNuclideBase, nuclideBases.DummyNuclideBase)
    ):
        return mapper.expand(armiName)  # warns, and leaves it out
    if isinstance(base, nuclideBases.NaturalNuclideBase):
        return mapper.expand(armiName)
    symbol = openmc.data.ATOMIC_SYMBOL[base.z]
    return [(f"{symbol}{base.a}" + (f"_m{base.state}" if base.state else ""), 1.0)]


def toArmiName(gnds):
    """The ARMI name of an OpenMC nuclide, or None if ARMI does not know it."""
    z, a, m = openmc.data.zam(gnds)
    base = nuclideBases.byAAAZZZSId.get(f"{a}{z:03d}{m}")
    return base.name if base is not None else None


def writeBack(burnable, densitiesById):
    """Set each depleted component's number densities and its block's burnup."""
    dropped = {}
    blocks = set()
    for mat, c, _volume in burnable:
        densities = {}
        for nuc, density in densitiesById[str(mat.id)].items():
            name = toArmiName(nuc)
            if name is None:
                dropped[nuc] = dropped.get(nuc, 0.0) + density
                continue
            densities[name] = densities.get(name, 0.0) + density
        c.setNumberDensities(densities)
        blocks.add(c.parent)
    for b in blocks:
        if b.p.molesHmBOL:
            b.p.percentBu = 100.0 * (1.0 - b.getHMMoles() / b.p.molesHmBOL)
    if dropped:
        worst = sorted(dropped, key=lambda nuc: dropped[nuc], reverse=True)[:5]
        runLog.warning(
            f"OpenMC depletion: {len(dropped)} nuclides ARMI does not know were left out of the "
            f"components, most abundant {worst}.",
            single=True,
        )
