Installation
============
Prerequisites:

- `CMake <http://www.cmake.org>`_ 3.4+
- A `C++11 compiler <http://en.cppreference.com/w/cpp/compiler_support>`_ (preferably GCC 4.8+ or Clang 3.3+)
- The `Boost <http://www.boost.org>`_ C++ libraries 1.50+, including runtime components  ``filesystem`` and ``program_options``
- (optional) The `HDF5 <http://www.hdfgroup.org/HDF5>`_ C++ library

All these dependencies are readily available on any operating system.
Some example installation commands for a few Linux distributions:

Ubuntu::

   apt-get install cmake g++ libboost-dev libboost-{filesystem,program-options}-dev libhdf5-dev

Fedora::

   yum install cmake gcc-c++ boost boost-{filesystem,program_options} hdf5 hd5-devel

Arch::

   pacman -S cmake gcc boost boost-libs hdf5

After installing the dependencies, download the `latest release <https://github.com/Duke-QCD/trento/releases/latest>`_ or `clone the repository <https://github.com/Duke-QCD/trento>`_, then compile and install T\ :sub:`R`\ ENTo through the standard CMake sequence::

   mkdir build && cd build
   cmake ..
   make install

This will install the compiled binary to ``~/.local/bin/trento``.
If you do not want this to happen, run ``make`` instead of ``make install`` and the binary will be left at ``build/src/trento``.
The remainder of this document assumes ``trento`` is in your ``PATH``.

The code is `continuously tested <https://travis-ci.org/Duke-QCD/trento>`_ on Ubuntu with GCC and Clang.
It should run just as well on any Linux distribution or OS X, and probably on Windows.
Other compilers should work but may require modifying the compiler flags.
