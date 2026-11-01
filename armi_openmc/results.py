"""
Reading a core statepoint and putting the results on the ARMI reactor.

OpenMC tallies are per source neutron. They become absolute rates once scaled by the source
rate that makes the tallied fission energy equal the reactor power::

    S [n/s] = P [W] / (sum of kappa-fission over the model [eV/n] * e [J/eV])

For each block the tallied totals are then divided by the number of copies of the block in the
model (an unfolded partial core holds each assembly several times) and by the block's symmetry
factor (ARMI stores only its share of a block cut by a symmetry line), which leaves exactly the
part of the block ARMI holds.
"""

import re
from dataclasses import dataclass

import numpy as np
import openmc
from armi import runLog

from armi_openmc.energy import toArmiOrder

EV_TO_J = 1.602176634e-19


@dataclass
class Tallied:
    mean: np.ndarray
    std: np.ndarray


class CoreResults:
    """The results of one core run, read from its statepoint."""

    def __init__(self, coreModel, statepointPath):
        self.coreModel = coreModel
        self.statepointPath = statepointPath
        with openmc.StatePoint(statepointPath) as sp:
            self.keff = sp.keff.nominal_value
            self.keffStd = sp.keff.std_dev
            self.totals = self._read(sp, "block totals")
            self.spectrum = self._read(sp, "block spectrum")
            self.pinTallies = {
                t.name: t.mean.ravel() for t in sp.tallies.values() if t.name.startswith("pins ")
            }
            self.entropy = np.array(sp.entropy) if sp.entropy is not None else None

    @staticmethod
    def _read(sp, name):
        tally = sp.get_tally(name=name)
        return Tallied(tally.mean, tally.std_dev)

    def sourceRate(self, power):
        """Neutrons per second that produce ``power`` watts in the modeled domain."""
        kappa = self.totals.mean[:, 0, 1].sum()
        if kappa <= 0.0:
            raise ValueError(
                "The OpenMC model tallied no fission energy; cannot normalize to power."
            )
        return power / (kappa * EV_TO_J)

    def apply(self, r, totalPower):
        """
        Set keff, block power, power density and flux on the reactor.

        Parameters
        ----------
        r : Reactor
        totalPower : float
            Full-core thermal power in W.
        """
        cm = self.coreModel
        core = r.core
        core.p.keff = self.keff
        core.p.openmcKeffStdDev = self.keffStd
        source = self.sourceRate(totalPower * cm.domainFraction)
        ng = len(cm.energyBounds) - 1
        cells = list(cm.blockCells)
        flux = self.totals.mean[:, 0, 0]
        kappa = self.totals.mean[:, 0, 1]
        mg = self.spectrum.mean[:, 0, 0].reshape(len(cells), ng)
        for idx, b in enumerate(cells):
            share = 1.0 / cm.copies.get(b, 1)
            if not cm.cut:
                share /= b.getSymmetryFactor()
            volume = b.getVolume()
            b.p.power = kappa[idx] * source * EV_TO_J * share
            b.p.pdens = b.p.power / volume
            b.p.mgFlux = toArmiOrder(mg[idx]) * source * share
            b.p.flux = flux[idx] * source * share / volume
            b.p.fluxPeak = b.p.flux
        if self.pinTallies:
            self._applyPinPowers(source, cells)
        runLog.info(
            f"OpenMC keff = {self.keff:.5f} +/- {self.keffStd:.5f}; "
            f"core power {sum(b.p.power for b in cells):.4e} W in {len(cells)} blocks"
        )

    def _applyPinPowers(self, source, cells):
        """Pin powers onto ``linPowByPin`` (W/cm), ordered like ``Block.getPinLocations``."""
        cm = self.coreModel
        byBlockCell = {cell.id: b for b, cell in cm.blockCells.items()}
        pinPower = {b: {} for b in cells}
        geometry_ = cm.model.geometry
        geometry_.determine_paths()
        latticeOf = {b: m.pinLattice for b, m in cm.blockModels.items() if m.pinLattice is not None}
        for cell in geometry_.get_all_cells().values():
            values = self.pinTallies.get(f"pins {cell.id}")
            if values is None:
                continue
            for path, value in zip(cell.paths, values, strict=True):
                b = next(
                    (byBlockCell[int(c)] for c in _CELLS.findall(path) if int(c) in byBlockCell),
                    None,
                )
                if b is None or b not in latticeOf:
                    continue
                lattice = latticeOf[b]
                found = [m for m in _LATTICE.finditer(path) if int(m.group(1)) == lattice.id]
                if not found:
                    continue
                index = tuple(int(v) for v in found[-1].group(2).split(","))
                ij = cm.blockModels[b].pinIndices.get(lattice.get_universe_index(index))
                if ij is not None:
                    pinPower[b][ij] = pinPower[b].get(ij, 0.0) + value * source * EV_TO_J
        for b, powers in pinPower.items():
            if not powers:
                continue
            share = 1.0 / cm.copies.get(b, 1) / (1.0 if cm.cut else b.getSymmetryFactor())
            height = b.getHeight()
            locations = b.getPinLocations() if hasattr(b, "getPinLocations") else []
            ordered = []
            for loc in locations:
                i, j, _ = loc.getCompleteIndices()
                ordered.append(powers.get((i, j), 0.0) * share / height)
            if ordered:
                b.p.linPowByPin = np.array(ordered)
            b.p.openmcPinPowers = np.array([powers[k] * share for k in sorted(powers)])
            b.p.openmcPinIndices = np.array(sorted(powers))


_CELLS = re.compile(r"c(\d+)")
_LATTICE = re.compile(r"l(\d+)\(([-\d,]+)\)")
