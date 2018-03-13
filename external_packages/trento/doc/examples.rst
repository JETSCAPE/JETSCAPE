Examples
========
Run a thousand lead-lead events using default settings and save the event data to file::

   trento Pb Pb 1000 > PbPb.dat

Run proton-lead events with a larger cross section (for the higher beam energy) and also compress the output::

   trento p Pb 1000 --cross-section 7.1 | gzip > pPb.dat.gz

Suppress printing to stdout and save events to HDF5::

   trento p Pb 1000 --cross-section 7.1 --quiet --output events.hdf

Uranium-uranium events at RHIC (smaller cross section) using short options::

   trento U U 1000 -x 4.2

Deformed gold-gold with an explicit nucleon width::

   trento Au2 Au2 1000 -x 4.2 -w 0.6

Simple sorting and selection (e.g. by centrality) can be achieved by combining standard Unix tools.
For example, this sorts by centrality (multiplicity) and selects the top 10%::

   trento Pb Pb 1000 | sort -rgk 4 | head -n 100

Working with Python
-------------------
The `scientific Python stack <https://www.scipy.org>`_ is ideal for analyzing T\ :sub:`R`\ ENTo data.

One way to load event properties into Python is to save them to a text file and then read it with ``np.loadtxt``.
Here's a nice trick to avoid the temporary file:

.. code:: python

   import subprocess
   import numpy as np

   with subprocess.Popen('trento Pb Pb 1000'.split(), stdout=subprocess.PIPE) as proc:
      data = np.array([l.split() for l in proc.stdout], dtype=float)

Now the ``data`` array contains the event properties.
It can be sorted and selected using numpy indexing, for example to sort by centrality as before:

.. code:: python

   data_sorted = data[data[:, 3].argsort()[::-1]]
   central = data_sorted[:100]

Text files are easily read by ``np.loadtxt``.
The header will be ignored by default, so this is all it takes to read and plot a profile:

.. code:: python

   import matplotlib.pyplot as plt

   profile = np.loadtxt('events/0.dat')
   plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)

Reading HDF5 files requires `h5py <http://www.h5py.org>`_.
``h5py`` file objects have a dictionary-like interface where the keys are the event names (``event_0``, ``event_1``, ...) and the values are HDF5 datasets.
Datasets can implicitly or explicitly convert to numpy arrays, and the ``attrs`` object provides access to the event properties.
Simple example:

.. code:: python

   import h5py

   # open an HDF5 file for reading
   with h5py.File('events.hdf', 'r') as f:
      # get the first event from the file
      ev = f['event_0']

      # plot the profile
      plt.imshow(ev, interpolation='none', cmap=plt.cm.Blues)

      # extract the profile as a numpy array
      profile = np.array(ev)

      # read event properties
      mult = ev.attrs['mult']
      e2 = ev.attrs['e2']

      # sort by centrality
      sorted_events = sorted(f.values(), key=lambda x: x.attrs['mult'], reverse=True)
