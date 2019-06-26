====================
T\ :sub:`R`\ ENTo 3D
====================

*Three-dimensional Extended Reduced Thickness Event-by-event Nuclear Topology*

A three-dimensional initial condition model for high-energy nuclear collisions. 
For a documentation of the original boost-invariant T\ :sub:`R`\ ENTo model, please read the docs at `qcd.phy.duke.edu/trento <http://qcd.phy.duke.edu/trento>`_. For model details please refer to `PRC 92 011901 <https://doi.org/10.1103/PhysRevC.92.011901>`_, `PRC 96 044912 <https://doi.org/10.1103/PhysRevC.96.044912>`_.

|

   |image1| |image2|

   **Figure 1**: *An event from* T\ :sub:`R`\ ENTo 3D. *Left: the event projected on to* :math:`x-y` *plane at midrapidty; right: same event projected on to* *y*\ -η *plane.*

|
|

**Installation**

First, install the dependencies: a C++11 compiler, Boost, GSL and HDF5. Then download the source code and compile with CMake:

.. code-block:: shell

   mkdir build && cd build
   cmake ..
   make install

|
|

**Usage**: 

T\ :sub:`R`\ ENTo 3D keeps all the options in the original boost-invariant model with additional options listed in Table 1. *Note that to work in three-dimensional mode*, :code:`--eta-max` *must be set to a positive value.*

Examples:

* Genereate 10 Pb Pb events at *s*\ :sup:`1/2` = 5020 GeV in boost-invariant mode and output to a :code:`.hdf5` file

.. code-block:: shell

   trento3d Pb Pb 10 -e 5020 -o PbPb.hdf5

* Same as above but generate three-dimensional initial condition with -10<η<10, *d*\ η=0.2.

.. code-block:: shell

   trento3d Pb Pb 10 -e 5020 --eta-max=10.0 --eta-step=0.2 -o PbPb.hdf5

* Use the absolute-skewness parametrization (see Table 2) instead of the relative-skewness parametrization, with skew coefficient γ\ :sub:`0`\ = 1.0.

.. code-block:: shell

   trento3d Pb Pb 10 -e 5020 -r 2 -t 1.0 --eta-max=10.0 --eta-step=0.2 -o PbPb.hdf5

.. csv-table:: **Table 1**: Additional program options
   :header: "Options", "Default", "Description"
   :widths: 10, 10, 35
   :align: center

   "-m, --mean-coeff", 1.0 (float>0), "rapidity mean coefficient *μ*\ :sub:`0` "
   "-s, --std-coeff", 3.0 (float>0), "rapidity std coefficient *σ*\ :sub:`0`"
   "-t, --skew-coeff", 0.0 (float>0), "rapidity skew coefficient *γ*\ :sub:`0`"
   "-r, --skew-type", 1 (int), "
					1 = relative skewness

					2 = absolute skewness
			
					else = no skewness"
   "-j, --jacobian", 0.8 (float>0), "<\ *p*\ :sub:`t`\ /\ *m*\ :sub:`t`\ > used in Jacobian *dy/d*\ η"
   "-e, --beam-energy", 2760 (float>0), "collision beam energy  *s*\ :sup:`1/2` [GeV], initializes cross section"
   "--xy-max",  10.0 (float) , "transverse x [fm] and y [fm] maximum (x,y grid from -max to +max)"
   "--xy-step",  0.2 (float), "transverse x [fm] and y [fm] step size"
   "--eta-max",  0.0 (float) , "space-time rapidity maximum (η grid from -max to +max)"
   "--eta-step",  0.5 (float), "space-time rapidity step size"

|
|

**Longitudinal extension**: 

T\ :sub:`R`\ ENTo 3D reproduces T\ :sub:`R`\ ENTo at midrapidity (η=0) exactly. At finite space-time rapidity, the entropy production is the product of its midrapidity value and a longitudinal profile function that varies at each transvese location. The profile is characterized by its first η-cumulants: mean, stadard deviationa and skewness. They are parametrized in terms of nuclear thickness function:

.. csv-table:: **Table 2**: cumulant parametrization
   :header: "Cumulants", "Parametrization"
   :widths: 15, 30
   :align: center

   "mean", "\ *μ*\ :sub:`0`\ /2 log [(\ *T*\ :sub:`A` *e*\ :sup:`Y`\ + \ *T*\ :sub:`B` *e*\ :sup:`-Y` ) / (\ *T*\ :sub:`A` *e*\ :sup:`-Y`\ + \ *T*\ :sub:`B` *e*\ :sup:`Y`\ )], *Y* is the beam rapidity"
   "standard deviation", *σ*\ :sub:`0`
   "skewness",  "Relative skewness, *γ*\ :sub:`0` (\ *T*\ :sub:`A`\ - \ *T*\ :sub:`B`\  )/(\ *T*\ :sub:`A`\ + \ *T*\ :sub:`B`\ )
   
   Absolute skewness, *γ*\ :sub:`0` (\ *T*\ :sub:`A`\ - \ *T*\ :sub:`B`\  )"

.. |image1| image:: doc/_static/event.png
   :width: 30%
.. |image2| image:: doc/_static/event-eta.png
   :width: 30%

