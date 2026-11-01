"""
Checks of a generated cross section library against the OpenMC run that made it.

With flux-volume weighted constants an infinite lattice keeps its neutron balance exactly, so the
eigenvalue of the homogenized infinite medium built from the ISOTXS has to reproduce the OpenMC
k-infinity of the block to within statistics. That tests the whole reduction: densities, flux
normalization, scattering matrix orientation and band storage, fission spectra, energy order.
"""

import numpy as np
from armi.nuclearDataIO import xsCollections


def macroscopic(lib, block):
    """ARMI's macroscopic cross sections of ``block`` from ``lib``."""
    creator = xsCollections.MacroscopicCrossSectionCreator()
    return creator.createMacrosFromMicros(lib, block)


def fluxWeightedChi(lib, block, flux):
    """The block's fission spectrum with every nuclide weighted by its actual fission source."""
    numerator = np.zeros(lib.numGroups)
    densities = block.getNumberDensities()
    for nuc in lib.getNuclides(block.getMicroSuffix()):
        m = nuc.micros
        source = densities.get(nuc.name, 0.0) * float(
            np.dot(m.neutronsPerFission * m.fission, flux)
        )
        numerator += m.chi * source
    return numerator / numerator.sum() if numerator.sum() > 0 else numerator


def infiniteMediumK(macros, chi=None):
    """
    k-infinity of a homogeneous medium from group constants.

    Solves ``(diag(total) - S) phi = (1/k) chi (nuSigF . phi)`` with ``S[to, from]`` the
    production-weighted scattering matrix; ``chi`` is the macroscopic fission-source spectrum.
    """
    total = np.ravel(macros.total[:, 0]) if np.ndim(macros.total) == 2 else np.asarray(macros.total)
    scatter = macros.totalScatter
    scatter = scatter.toarray() if hasattr(scatter, "toarray") else np.asarray(scatter)
    nuSigF = np.asarray(macros.nuSigF)
    chi = np.asarray(macros.chi if chi is None else chi)
    if chi.sum() <= 0:
        return 0.0
    chi = chi / chi.sum()
    loss = np.diag(total) - scatter
    fissionMatrix = np.outer(chi, nuSigF)
    eig = np.linalg.eigvals(np.linalg.solve(loss, fissionMatrix))
    return float(np.max(eig.real))
