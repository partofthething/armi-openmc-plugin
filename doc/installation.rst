Installation
============

The plugin needs three things: ARMI (the partofthething fork), the OpenMC Python API, and an
OpenMC solver with nuclear data.

Getting the code
----------------

.. code-block:: bash

   git clone https://github.com/partofthething/armi-openmc-plugin.git
   cd armi-openmc-plugin

Python packages
---------------

PyPI has no OpenMC wheels. The package therefore installs the OpenMC Python API from its release
tag; it is pure Python, enough to build models and read statepoints. ``pyproject.toml`` pins ARMI
to the fork's ``main`` and OpenMC to ``v0.15.3`` through ``[tool.uv.sources]``:

.. code-block:: bash

   uv venv
   uv pip install -e ".[test]"

To develop against a local ARMI checkout, install it editable *afterwards*, because installing
this package pulls ARMI from git and replaces whatever ARMI was there:

.. code-block:: bash

   uv pip install --no-deps -e ../armi

The solver
----------

``openmcRunner: auto`` uses an ``openmc`` executable on the ``PATH`` and falls back to docker.
The container route needs nothing installed but docker:

.. code-block:: bash

   docker pull openmc/openmc:latest

The plugin mounts the run directory and the directory holding the nuclear data at the same paths
inside the container as outside, so the paths written into ``model.xml`` resolve either way. It
runs the container as your user, so outputs belong to you.

Nuclear data
------------

Any OpenMC HDF5 library works; the benchmarks here use the official ENDF/B-VIII.1 library from
https://openmc.org/official-data-libraries/. Point the ``openmcCrossSections`` setting, or the
``OPENMC_CROSS_SECTIONS`` environment variable, at its ``cross_sections.xml``.

Building these docs
-------------------

.. code-block:: bash

   uv pip install sphinx sphinx-book-theme
   python doc/make_figures.py          # refresh the benchmark figures from the example runs
   sphinx-build -b html doc doc/_build/html
