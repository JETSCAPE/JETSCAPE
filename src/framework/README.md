Test Skeleton of (potential) JetScape Framework

With cmake:
- on shell type: cmake .
- on shell type: make
and the you should have to test binaries in ./bin so test:
- ./bin/brickTest

If cmake found other libraries HepMC, ROOT or Pythia8, you might have to add the librry path's in the setup.csh script.
For sure works on Mac Os X 10.11.6.

---> Obsolete
To compile and use (Makefiles):
Create a lib and bin directory first (the Makefiles and the setup.csh script for the shared library assume that structure)
Create Library (statitc or shared):

cd src
make (This makes the shared libarary)
make -f Makefile.static (This makes a static Library)

Remark: Once we agree on if we want to use Makefiles or cmake (a simple example conf for cmake can be found in the CMakeConf directory)
the Makefiles will be extend and merged.

Go back to main directory and complie the test programm (./src/testJetScape.cc):

make (creates the test program using the shared library in ./lib/libJetScape.so)

in order to use it set the paths (using a csh: source setup.csh)

and execute

bin/testJetScape

If you created a static library and want to use it for the test program:

make static

bin/testJetScape

Remark: The Makefiles should detect if you have a linux or mac and should change the flags as necessary.
So far everything works for Mac OSX 10.11.6 and changes to Makefile on linux might be necessary 

To make a class documentation using doxygen:

doxygen JetScapeDoxy.conf
