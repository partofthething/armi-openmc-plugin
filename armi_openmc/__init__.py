"""
ARMI plugin for the OpenMC Monte Carlo code.

* ``neutronicsKernel: OpenMC``: whole-core eigenvalue runs, pin-detailed or homogenized, in
  hexagonal and Cartesian cores, with keff, block power and multigroup flux mapped back.
* ``xsKernel: OpenMC``: microscopic multigroup cross sections for each cross section group from
  a block or colorset model, written as ISOTXS.

Register :py:class:`armi_openmc.plugin.OpenMCPlugin` with an ARMI app, or run
``python -m armi_openmc``, whose app already has it.
"""

__version__ = "0.1.0"
