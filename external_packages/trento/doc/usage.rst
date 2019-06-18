Usage
=====
T\ :sub:`R`\ ENTo has a standard command-line interface.
The basic syntax is ::

   trento [options] projectile projectile [number-events = 1]

where the only required arguments are the two projectile names.
For example, ``trento Pb Pb 10`` would run ten lead-lead events.

The remaining optional arguments may be given in any order, before or after the projectiles.
Run ``trento --help`` for a brief summary of the options and see below for more detailed descriptions.

Specifying projectiles
----------------------
The ``projectile`` arguments take species abbreviations, e.g. ``p``, ``Pb``, etc.
The known species are

=========  ========  ============  ========
Symbol     Name      No. nucleons  Deformed
=========  ========  ============  ========
p          proton    1             ---
d          deuteron  2             ---
Cu         copper    63            no
Cu2        copper    63            yes
Xe         xenon     129           no
Au         gold      197           no
Au2        gold      197           yes
Pb         lead      208           no
U, U2, U3  uranium   238           yes
=========  ========  ============  ========

For the deuteron, nucleon positions are sampled from the Hulthén wavefunction;
for the heavy nuclei, positions are sampled from a `Woods-Saxon <https://en.wikipedia.org/wiki/Woods%E2%80%93Saxon_potential>`_ distribution, either spherically symmetric or deformed as indicated.
Copper and gold are slightly deformed—slightly enough that a symmetric distribution is a reasonable approximation—therefore both symmetric (``Cu``, ``Au``) and deformed (``Cu2``, ``Au2``) versions are provided, where both versions have the same nuclear radius and surface thickness.
There is no consensus on the uranium Woods-Saxon parameters, so three commonly used sets are provided:

======  ====  ====  =====  =====
Symbol  *R*   *a*   |b2|   |b4|
======  ====  ====  =====  =====
U       6.81  0.60  0.280  0.093
U2      6.86  0.42  0.265  0
U3      6.67  0.44  0.280  0.093
======  ====  ====  =====  =====

.. |b2| replace:: *β*\ :sub:`2`
.. |b4| replace:: *β*\ :sub:`4`

The ``U`` and ``U2`` sets are given in this recent `overview of particle production from PHENIX <http://inspirehep.net/record/1394433>`_.
All other Woods-Saxon parameters (including ``U3``) and the Hulthén wavefunction parameters are from the `PHOBOS Glauber model <http://inspirehep.net/record/1310629>`_.

For Woods-Saxon nuclei, ``trento`` can impose a minimum nucleon-nucleon distance.
See the :ref:`nucleon-min-dist <nucleon-min-dist>` option.

In addition, ``trento`` can read :ref:`arb-configs` saved in HDF5 files.

General options
---------------
These are general options that don't fit in any other category.

``-h, --help``
   Show the help message and exit.

``--version``
   Print version number and exit.

``--bibtex``
   Print bibtex entry and exit.

``-c, --config-file FILE``
   Path to configuration file (see :ref:`config-files` below).
   May be given multiple times.


Output options
--------------
The default output mode is to print event-by-event properties to stdout, in the following order::

   event_number impact_param npart mult e2 e3 e4 e5

with one line for each event, where

- ``event_number`` is an integer counter,
- ``impact_param`` is the collision impact parameter,
- ``npart`` is the number of nucleon participants,
- ``mult`` is the total initial entropy, and
- the ``en`` are the eccentricity harmonics ɛ\ :sub:`n`.

This format is designed for easy parsing, redirection to files, etc.
The output may be disabled with the ``-q/--quiet`` option.

By default, the actual initial entropy profiles (grids) are not output.
There are two available output formats: text and HDF5 (if compiled).

In text mode, each event is written to a separate text file as a standard block-style grid, along with a commented header containing the event properties, like this::

   # event 0
   # b     = 2.964077155
   # npart = 380
   # mult  = 168.603282
   # e2    = 0.01953253866
   # e3    = 0.08961920965
   # e4    = 0.1101683349
   # e5    = 0.1727159106

The header may be disabled with the ``--no-header`` option.

HDF5 is a high-performance, cross-platform binary format for large numerical datasets.
Libraries are available in `most languages <https://en.wikipedia.org/wiki/Hierarchical_Data_Format#Interfaces>`_.
HDF5 is significantly faster than text output:
writing an event to a text file usually takes much longer than computing the actual event;
writing to HDF5 incurs only a small overhead.
Therefore, HDF5 is the recommended output format.

In HDF5 mode, all events are written to a single file with each event in a separate HDF5 dataset.
Event properties are written to each dataset as HDF5 attributes with names ``b``, ``npart``, ``mult``, ``e2``, etc.

``-q, --quiet``
   Disable printing event properties to stdout.
   Since both text and HDF5 output contain the event properties, it's often desirable to specify this option along with the output option.

``-o, --output PATH``
   Path to output events.
   If the path has an HDF5-like extension (``.hdf5``, ``.hdf``, ``.hd5``, ``.h5``), then all events will be written to that HDF5 file.
   Otherwise, the path is interpreted as a directory and events will be written to numbered text files in the directory.

   For text output, the directory will be created if it does not exist.
   If it does already exist, it must be empty (this is to avoid accidentally overwriting files or spewing thousands of files into an already-used location).

   For HDF5 output, the file must not already exist.
   Each event will be written as a numbered dataset in the file, and the standard event properties will be written as dataset attributes.

   Example:

   - ``--output events`` will write to text files ``events/0.dat``, ``events/1.dat``, ...
   - ``--output events.hdf`` will write to HDF5 file ``events.hdf`` with dataset names ``event_0``, ``event_1``, ...

``--no-header``
   Disable writing event headers to text files.

Physical options
----------------
These options control the physical behavior of the model.

.. warning::

   The physical options have reasonable defaults, however **the defaults are not in any way a best-fit to experimental data**.
   They are simply round numbers.
   It is entirely expected that the ideal parameters will change depending on the beam energy.
   In particular, **the cross section must be explicitly set for each beam energy**.

``-p, --reduced-thickness FLOAT``
   Reduced thickness parameter *p*.
   The reduced thickness is defined as the `generalized mean <https://en.wikipedia.org/wiki/Generalized_mean>`_ of participant nuclear thickness

   .. math::

      T_R(p; T_A, T_B) = \biggl( \frac{T_A^p + T_B^p}{2} \biggr)^{1/p}

   The default is *p* = 0, which corresponds to the geometric mean.

``-k, --fluctuation FLOAT``
   `Gamma distribution <https://en.wikipedia.org/wiki/Gamma_distribution>`_ shape parameter *k* for nucleon fluctuations.
   Fluctuations are sampled from a gamma distribution with the scale parameter fixed so that the mean is one:

   .. math::

      P_k(x) = \frac{k^k}{\Gamma(k)} x^{k-1} e^{-kx}

   The default is *k* = 1, which corresponds to an exponential distribution.
   For small *k*, the distribution has a long tail, leading to large fluctuations.
   For large *k*, the distribution becomes a narrow Gaussian, and eventually a delta function for very large values.

``-w, --nucleon-width FLOAT``
   Gaussian nucleon width in fm:

   .. math::

      T_\text{nucleon}(x, y) = \frac{1}{2\pi w^2} \exp\biggl( -\frac{x^2 + y^2}{2w^2} \biggr)

   The default is 0.5 fm.

.. _nucleon-min-dist:

``-d, --nucleon-min-dist FLOAT``
   Minimum nucleon-nucleon distance (fm) for Woods-Saxon nuclei (spherical and deformed).
   When nonzero, if a sampled nucleon lands too close to a previously sampled nucleon, its angular position is resampled until it lands far enough away.
   The radius is *not* resampled, since this would effectively modify the Woods-Saxon distribution.

   If a nucleon cannot be placed after a reasonable number of retries, the algorithm gives up and leaves the nucleon at the last sampled position.
   The failure rate is negligible for minimum distances of ~1 fm and below;
   it reaches roughly 1% at 1.7 fm for spherical nuclei and 1.5 fm for deformed.

   The default is zero (no minimum distance).

   .. versionadded:: 1.4

``-x, --cross-section FLOAT``
   Inelastic nucleon-nucleon cross section |snn| in |fm2|.
   The default is 6.4 fm\ :sup:`2`, the approximate experimental value at LHC Pb+Pb energy, √s = 2.76 TeV.
   Here are some measurements of the cross section at common beam energies (all have approximately 0.5 |fm2| uncertainty):

   +---------+---------------+---------------+
   |√s [TeV] | |snn| [|fm2|] | ref.          |
   +=========+===============+===============+
   |0.200    | 4.23          | `1509.06727`_ |
   +---------+---------------+---------------+
   |         | 6.4           | `1108.6027`_  |
   + 2.76    +---------------+---------------+
   |         | 6.28          | `1208.4968`_  |
   +---------+---------------+---------------+
   |5.02     | 7.0           | `1210.3615`_  |
   +---------+---------------+---------------+
   |7        | 7.32          | `1208.4968`_  |
   +---------+---------------+---------------+

.. |snn| replace:: σ\ :sub:`NN`
.. |fm2| replace:: fm\ :sup:`2`
.. _1108.6027: https://inspirehep.net/record/925723
.. _1210.3615: https://inspirehep.net/record/1190545
.. _1208.4968: https://inspirehep.net/record/1181770
.. _1509.06727: https://inspirehep.net/record/1394433

``-n, --normalization FLOAT``
   Overall normalization factor.
   The default is 1.

``--b-min FLOAT``
   Minimum impact parameter.
   The default is zero.

``--b-max FLOAT``
   Maximum impact parameter.
   The default is to run minimum-bias collisions for the given collision system.

   To run at fixed impact parameter, give the same value for both the min and the max.

``--random-seed POSITIVE_INT``
   Primarily for testing and debugging.

Grid options
------------
The thickness functions are discretized onto a square *N* × *N* grid centered at (0, 0).
The grid can have a dramatic effect on code speed and precision, so should be set carefully.
Computation time is roughly proportional to the number of grid cells (i.e. *N*\ :sup:`2`).

``--grid-max FLOAT``
   *x* and *y* maximum of the grid in fm, i.e. the grid extends from -max to +max.
   The default is 10 fm, large enough to accommodate all collision systems.
   However, this should be set as small as possible, since an unnecessarily large grid slows down the code.
   For anything but uranium-uranium, 9 fm is sufficient.
   For pp and pA, 3 fm is usually a good choice.

``--grid-step FLOAT``
   Size of grid cell in fm.
   The default is 0.2 fm, sufficient to achieve ~99.9% precision for the event properties.
   This can reasonably be increased as far as the nucleon width; beyond that and precision suffers significantly.

The grid will always be a square *N* × *N* array, with *N* = ceil(2*max/step).
So e.g. the default settings (max = 10 fm, step = 0.2 fm) imply a 100 × 100 grid.
The ceiling function ensures that the number of steps is always rounded up, so e.g. given max = 10 fm and step 0.3 fm, the grid will be 67 × 67.
In this case, the actual grid max will be marginally increased (max = nsteps*step/2).

Regardless of the collision system, the code will always approximately center the overlap region on the grid.

.. _config-files:

Configuration files
-------------------
.. highlight:: ini

All options may be saved in configuration files and passed to the program via the ``-c, --config-file`` option.
Config files follow a simple ``key = value`` syntax, and lines beginning with a ``#`` are comments.
The key for each option is its long option without the ``--`` prefix.
Here's an example including all options::

   # specify the projectile option twice
   projectile = Pb
   projectile = Pb
   number-events = 1000

   # don't print event properties to stdout, save to HDF5
   quiet = true
   output = PbPb.hdf

   reduced-thickness = 0
   fluctuation = 1
   nucleon-width = 0.5
   cross-section = 6.4
   normalization = 1

   # leave commented out for min-bias
   # b-min =
   # b-max =

   grid-max = 10
   grid-step = 0.2

Multiple config files can be given and they will be merged, so options can be separated into modular groups.
For example, one could have a file ``common.conf`` containing settings for all collision systems and files ``PbPb.conf`` and ``pp.conf`` for specific collision systems::

   # common.conf
   reduced-thickness = 0.2
   fluctuation = 1.5
   nucleon-width = 0.6

   # PbPb.conf
   projectile = Pb
   projectile = Pb
   number-events = 10000
   grid-max = 9

   # pp.conf
   projectile = p
   projectile = p
   number-events = 100000
   grid-max = 3

.. highlight:: none

To be used like so::

   trento -c common.conf -c PbPb.conf
   trento -c common.conf -c pp.conf

If an option is specified in a config file and on the command line, the command line overrides.

.. _arb-configs:

Arbitrary nuclear configurations
--------------------------------
.. versionadded:: 1.3

``trento`` can read pre-generated nuclear configurations from HDF5 files.

The following files were created from publicly available data and can be input directly to ``trento``.
They are redistributed with permission from the authors.

- |3He| configurations are from the `PHOBOS Glauber model <https://tglaubermc.hepforge.org>`_, created by Joe Carlson at LANL (`ref <http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.70.743>`_).
- |197Au| and |208Pb| configurations including realistic nucleon-nucleon correlations were created by Massimiliano Alvioli (`ref 1 <http://inspirehep.net/record/820666>`_, `ref 2 <http://inspirehep.net/record/1082705>`_) and are available on `his website <http://users.phys.psu.edu/~malvioli/eventgenerator>`_.

If you use these configurations in your research, please cite the original authors.

=======  ===============  ===========  =======  ============================================
Species  File             No. configs  Size     sha1sum
=======  ===============  ===========  =======  ============================================
|3He|    He3.hdf_          13,699      484 KiB  ``a50c22ad8999db185e50fa513adf8100c29fba8c``
|197Au|  Au197.hdf_         1,820      4.2 MiB  ``9124eeab163bb2fbc6a919cb96efd44b99cac6be``
|208Pb|  Pb208_10k.hdf_    10,000       24 MiB  ``4d5c76cb4b5535538b57864a1287a4695abc29d1``
|208Pb|  Pb208_100k.hdf_  100,000      239 MiB  ``d67f7aca2b14f8c705a4bfa0a8aeedcd3a816f6e``
=======  ===============  ===========  =======  ============================================

.. |3He| replace:: :sup:`3`\ He
.. |197Au| replace:: :sup:`197`\ Au
.. |208Pb| replace:: :sup:`208`\ Pb

.. _He3.hdf: nuclear-configs/He3.hdf
.. _Au197.hdf: nuclear-configs/Au197.hdf
.. _Pb208_10k.hdf: nuclear-configs/Pb208_10k.hdf
.. _Pb208_100k.hdf: nuclear-configs/Pb208_100k.hdf

Pb208_10k.hdf_ contains the same data as the first 10,000 configurations in Pb208_100k.hdf_.
The smaller file is provided for convenience.

To use pre-generated configurations, specify a path to an appropriate file on the command line in place of a species abbreviation::

   trento path/to/file1.hdf path/to/file2.hdf

Filenames must have an HDF5-like extension (``.hdf5``, ``.hdf``, ``.hd5``, ``.h5``).
The files may be the same or different and may be mixed with standard species abbreviations.
For each event, ``trento`` will choose a random configuration from the file and apply a random three-dimensional rotation.
Hence, it is safe to run several events per pre-generated configuration.

For example, to run |3He|\ +Au events at RHIC, download He3.hdf_ and execute ::

   trento --cross-section 4.2 He3.hdf Au2

Remember to set the appropriate cross section for the desired beam energy!

To run custom configurations, make an HDF5 file containing a single dataset of shape ``(number_configs, number_nucleons, 3)``, where the first dimension corresponds to each configuration, the second dimension to each nucleon, and the third dimension to the (x, y, z) coordinates of each nucleon.
Note that ``trento`` will read the file as single-precision floats, not doubles.

.. highlight:: python

The easiest way to write an HDF5 file is with `h5py <http://www.h5py.org>`_::

   import numpy as np
   import h5py

   # generate random data for 10 configs of a nucleus with 100 nucleons
   configs = np.random.uniform(-1, 1, (10, 100, 3))

   with h5py.File('nuclear_configs.hdf') as f:
      # the name of the dataset does not matter as long as there is only one
      f.create_dataset('configs', data=configs, dtype=np.float32)

