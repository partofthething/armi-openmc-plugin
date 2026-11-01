"""ARMI energy group structures as OpenMC energy filter bounds."""

import numpy as np
from armi.physics.neutronics import energyGroups

#: Lower edge of the last group. ARMI structures list only upper bounds; OpenMC data starts at
#: 1e-5 eV.
MINIMUM_ENERGY = 1.0e-5


def upperBounds(name):
    """Group upper bounds in eV, highest group first (ARMI order)."""
    return np.array(energyGroups.getGroupStructure(name), dtype=float)


def energyFilterBounds(name):
    """Ascending bounds for an :py:class:`openmc.EnergyFilter`, ``numGroups + 1`` long."""
    upper = upperBounds(name)
    return np.concatenate(([MINIMUM_ENERGY], upper[::-1]))


def toArmiOrder(values, axis=-1):
    """Flip an array indexed by ascending OpenMC energy bin into ARMI's fast-first order."""
    return np.flip(np.asarray(values), axis=axis)
