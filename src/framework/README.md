Test Skeleton of (potential) JetScape Framework

With cmake:

- on shell type: cmake .
- on shell type: make

and the you should have be able to use the test binaries in ./bin.
So, you can start the testing now:

- ./bin/brickTest

This "simulates" our brick test. The temperature from the brick is available in the jet energy loss modules. The jet energyloss module consists of two "test" modules, Matter performing a random in time democratic split and Martini is doing nothing. Furthermore in "Matter" I added random "new graph roots" for testing (simulating scattering with medium partons). The switching criteria is arbitrarily set to 5GeV in pt for testing. An ascii output file is created which you can read in with the provided reader:

- ./readerTest

which reads in the generated showers does some DFS search and shows the output. You can generate an output graph format which can be easily visualized using graphViz or other tools like Gephi (GUI for free for Mac) or more adanvanced, graph-tools (Python) and many more. Furthermore, as a "closure" test, I used the FastJet core packackage (compiled in our JetScape library) to perform a simple jetfinding (on the "final" partons, in graph language, incoming partons in a vertex with no outgoing partons/childs), and since the "shower" is perfectly collinear the jet pT is identical to the hard process parton pT (modulo some random new partons/roots in the final state, see above).  

If cmake found other libraries HepMC, ROOT or Pythia8, you might have to add the library path's in the setup.csh script.
For sure works on Mac Os X 10.11.6.

To make a class documentation using doxygen:

doxygen JetScapeDoxy.conf

Getting started on a Mac:

- Install Xcode and command-line tools

For further packages needed (like CMake, Pythia, ROOT, GraphViz) I recommend homebrew for Mac:

- Install homebrew (after you install Xcode)

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

- Install CMake via homebrew type: 

brew install cmake

- Install doxygen:

brew install doxygen

- Tap some repos for further sceintific packages: 

brew tap davidchall/hep

brew tap homebrew/science

- Install Pythia8 for example:

brew install pythia8

- Install graphViz:

brew install graphviz --with-app

- Install root6

brew install root6

(and most of other packackes we need)

Remark: So far on brew HepMC is only available in version 2, version 3 is required for the current code, nicer wirter interfaces to root for example. So one has to install it from: http://hepmc.web.cern.ch/hepmc/

