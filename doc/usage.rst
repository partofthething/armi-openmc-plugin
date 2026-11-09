Using the plugin
================

The plugin registers with any ARMI application. ``python -m armi_openmc`` is ARMI's command
line with the plugin already loaded; an application of your own registers
``armi_openmc.plugin.OpenMCPlugin`` the usual way. Or let the case choose it, with the fork's
``plugins`` setting, and run plain ARMI:

.. code-block:: yaml

    settings:
      plugins:
        - armi_openmc.plugin.OpenMCPlugin
      neutronicsKernel: OpenMC

.. code-block:: bash

    python -m armi run case.yaml

Global flux
-----------

Set ``neutronicsKernel: OpenMC`` and run the case as usual:

.. code-block:: yaml

   settings:
     loadingFile: blueprints.yaml
     power: 400000000.0
     neutronicsKernel: OpenMC
     openmcCrossSections: /path/to/endfb-viii.1-hdf5/cross_sections.xml
     openmcParticles: 40000
     openmcBatches: 200
     openmcInactiveBatches: 50
     openmcBlockModel: heterogeneous
     openmcHeterogeneousBlockFlags: [fuel, control]   # empty: every block in detail

At every time node the interface builds the model, runs it in its own folder under
``openmcRunDir`` (``model.xml``, the statepoint, ``openmc.out`` and, with ``openmcPlotGeometry``,
PNG slices), and sets on the reactor:

* ``core.p.keff`` and ``core.p.openmcKeffStdDev``;
* block ``power``, ``pdens``, ``flux`` and ``mgFlux`` (volume-integrated, ARMI group order, in
  the ``groupStructure`` of the case or ``openmcGroupStructure``);
* with ``openmcPinPowers``, ``linPowByPin`` (W/cm, ordered like ``Block.getPinLocations``) and
  ``openmcPinPowers`` / ``openmcPinIndices`` (W per pin, with its block grid indices).

Results are normalized so that the fission energy tallied in the model equals the case's
``power`` (scaled to the modeled fraction of the core); see :doc:`modeling`.

OpenMC's batch output is echoed into the ARMI log as ``[openmc]`` lines, so a terminal, a log file
and the input editor's console all show the run converging.

Depletion
---------

Set ``openmcDepletion: true`` and the plugin burns the reactor with ``openmc.deplete``. At every
time node, before the flux solve, it depletes the reactor over the burn step that just ended:

1. it builds the whole-core model from the reactor as it stands, and marks depletable every
   material that is a single depletable component (``DEPLETABLE`` flag), with its volume;
2. it runs one ``openmc.deplete`` step (``deplete_step.py``, copied into
   ``<openmcRunDir>/deplete-c<cycle>n<node>``) with the step's full-power days and power, the
   chain (``openmcDepletionChain`` or ``OPENMC_CHAIN_FILE``) and the integrator
   (``openmcDepletionIntegrator``, CE/CM by default);
3. it writes the number densities at the end of the step back onto the components, every
   nuclide of the chain that ARMI knows (fission products and actinides alike), and updates each
   block's ``percentBu`` from its heavy metal.

The outage between cycles (``cycleLength * (1 - availabilityFactor)``) is decayed at zero power,
without transport. The inventory sent to OpenMC includes the nuclides the data library has no
cross sections for; OpenMC carries those through decay, and they come back onto ARMI, so nothing
is lost between steps.

``openmcDepletionReactionRates: flux`` (the default) tallies the flux in
``openmcDepletionGroups`` (SHEM-361) and collapses each nuclide's cross sections with it, with
fission and capture of the major actinides tallied directly; burned fuel carries hundreds of
nuclides, and tallying every reaction of each (``direct``) makes the transport several times
slower. On a 365-day step of a small LWR core the two agree within 0.4% on every core-total
inventory checked (U-235, Pu, Sm-149, Cs-137, Nd-148).

The depletion interface does not set keff or the power distribution: pair it with
``neutronicsKernel: OpenMC``, which then solves the depleted reactor at every node. Only
components that are their own OpenMC material can deplete, i.e. those of heterogeneous blocks.
Depletion needs OpenMC's Python API with its shared library: ``openmcRunner: auto`` uses
``openmcPython`` when it can ``import openmc.lib`` and the container otherwise.

.. code-block:: yaml

    neutronicsKernel: OpenMC
    openmcDepletion: true
    openmcDepletionChain: /data/chain_endfb71.xml
    cycles:
      - step days: [1, 4, 10, 15, 30, 60, 120, 240, 365, 365]
        availability factor: 1.0

A chain matched to a library can be made from its ENDF decay, fission yield and neutron
sublibraries with ``openmc.deplete.Chain.from_endf`` (the ``openmc/openmc`` image carries
ENDF/B-VII.1's in ``/root/endf-b-vii.1``).

Cross sections
--------------

Set ``xsKernel: OpenMC`` (with ``genXS: Neutron``) and the lattice physics interface makes the
ISOTXS for every cross section group the ARMI cross section group manager defines. A group whose
cross section settings name a ``driverID`` is modeled surrounded by its driver when
``openmcXSColorset`` is on. :doc:`crosssections` explains the models and the reduction.

Command line
------------

.. code-block:: bash

   python -m armi_openmc openmc-core case.yaml [--no-run] [--plot] [--label NAME]
   python -m armi_openmc openmc-xs   case.yaml [--xs-ids AA BA] [--compare OTHER_ISOTXS]
   python -m armi_openmc compare-xs  REFERENCE_ISOTXS OTHER_ISOTXS [--flux FLUX.json]

``openmc-core`` builds (and runs) the whole-core model of a case and writes ``results.json``
with keff, assembly and block powers, and pin powers. ``openmc-xs`` makes an ISOTXS outside a
full run. ``compare-xs`` compares two libraries nuclide by nuclide, collapsing with the group flux
the OpenMC library was made with (``<library>.flux.json``, written beside it).

The input editor
----------------

``python -m armi_openmc gui case.yaml`` opens ARMI's browser-based input editor with the plugin
loaded. The plugin's settings form their own **OpenMC** section, ordered, described, and shown
only when relevant (the executable only for the native runner, the block flags only for the
heterogeneous model).

.. figure:: _static/webui-settings.png
   :width: 100%
   :alt: The OpenMC settings section of the ARMI input editor

   The OpenMC settings of the 2D C5G7 case in the input editor.

**Run** (Ctrl+Enter) runs the open case with the same application, here with the plugin, and
follows its log in a console docked under the workspace, with the elapsed time and a Stop
button. Unsaved edits are saved first, since a run reads the files on disk.

.. figure:: _static/webui-run.png
   :width: 100%
   :alt: A C5G7 run in progress in the input editor's run console

   A C5G7 run started from the editor.
