Module reference
================
Each module encapsulates a unit of functionality in the parameter estimation project.
The following is a brief description of each module and its most important classes, functions, and variables.

.. note::

   Some downloaded data and results of long-running functions are saved on disk in the :file:`cache` directory, which will be created as needed.
   However, there is no logic for updating the cache if the code is modified, which means cache files may need to be occassionally manually deleted.

Package init file
-----------------
Source code: :ghlink:`src/__init__.py`

.. automodule:: src

.. autodata:: systems

.. autofunction:: parse_system

Experimental data
-----------------
Source code: :ghlink:`src/expt.py`

.. automodule:: src.expt

.. autoclass:: HEPData
   :members:

.. autofunction:: _data

.. autodata:: data
   :annotation: = <nested dict object>

.. autofunction:: cov

Model data
----------
Source code: :ghlink:`src/model.py`

.. automodule:: src.model

.. autoclass:: ModelData
   :members:

.. data:: data

   A nested dict of model data with the same structure as :data:`src.expt.data`.

Design
------
Source code: :ghlink:`src/design.py`

.. automodule:: src.design

.. autofunction:: generate_lhs

.. autoclass:: Design
   :members:

Emulator
--------
Source code: :ghlink:`src/emulator.py`

.. automodule:: src.emulator

.. autoclass:: Emulator
   :members:

MCMC
----
Source code: :ghlink:`src/mcmc.py`

.. automodule:: src.mcmc

.. autoclass:: Chain(path=Path('mcmc/chain.hdf'))
   :members:
   :exclude-members: map

.. autofunction:: credible_interval

Plots and figures
-----------------
Source code: :ghlink:`src/plots.py`

.. automodule:: src.plots
