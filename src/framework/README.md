Test Skeleton of (potential) JetScape Framework

With cmake:
- on shell type: cmake .
- on shell type: make
and the you should have be able to use the test binaries in ./bin.
Sonstart the testing:
- ./bin/brickTest
This "simulates" our brick test, so the temperature from the brick is available in the jet enerhyloss module. The jet energyloss module consists of two "test" modules, Matter performing a random in time democratic split and Martini is doung nothing. The switching criteria is arbitrarily set to 5GeV in pt for testing. An ascii output file is created which you can read in with the privided reader:
- ./readerTest
which reads in the generated showers does some DFS search and shows the output, as as can generate an output graph format which can be easily visualized using graphViz.

If cmake found other libraries HepMC, ROOT or Pythia8, you might have to add the library path's in the setup.csh script.
For sure works on Mac Os X 10.11.6.

To make a class documentation using doxygen:

doxygen JetScapeDoxy.conf
