# JETSCAPE1.0-RC1
The Jet Energy-loss Tomography with a Statistically and Computationally Advanced Program Envelope (JETSCAPE) collaboration


With cmake:

```bash
    mkdir build
    cd build
    cmake ..
```

the test binaries will be generated in the build folder
To run the tests, you need to stay at the framework directory and type :

```bash
    ./build/brickTest
```

This "simulates" our brick test. The temperature from the brick is available in the jet energy loss modules. The jet energyloss module consists of two "test" modules, Matter performing a random in time democratic split and Martini is doing nothing. Furthermore in "Matter" I added random "new graph roots" for testing (simulating scattering with medium partons). The switching criteria is arbitrarily set to 5GeV in pt for testing. An ascii output file is created which you can read in with the provided reader:

```bash
    ./build/readerTest
```

which reads in the generated showers does some DFS search and shows the output. You can generate an output graph format which can be easily visualized using graphViz or other tools like Gephi (GUI for free for Mac) or more adanvanced, graph-tools (Python) and many more. Furthermore, as a "closure" test, I used the FastJet core packackage (compiled in our JetScape library) to perform a simple jetfinding (on the "final" partons, in graph language, incoming partons in a vertex with no outgoing partons/childs), and since the "shower" is perfectly collinear the jet pT is identical to the hard process parton pT (modulo some random new partons/roots in the final state, see above).  

```bash
    ./build/hydroFromFileTest
```

This "simulates" a test by reading in hydrodynamic evolution from
an external file. The temperature and flow velocity at any given space-time
poisition is available in the jet energy loss modules.
The jet energyloss module is the same as the one in the brickTest.
The hydro file reader is from a 3rd-party program. Part of the code requires
the HDF5 library.

```bash
    ./build/MUSICTest
```

This "simulates" a test with MUSIC (a viscous hydrodynamic code).
The temperature and flow velocity at any given space-time position
is available in the jet energy loss modules.
The jet energyloss module is the same as the one in the brickTest.
MUSIC is treated as a 3rd party program outside the JETSCAPE framework.
It requires MPI and GSL libraries installed in the computer.

If cmake found other libraries HepMC, ROOT or Pythia8, you might have to add the library path's in the setup.csh script.
For sure works on Mac Os X 10.11.6.

To make a class documentation using doxygen:

```bash
doxygen JetScapeDoxy.conf
```

Getting started on a Mac:

- Install Xcode and command-line tools

For further packages needed (like CMake, Pythia, ROOT, GraphViz) I recommend homebrew for Mac:

- Install homebrew (after you install Xcode)

```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

- Install CMake via homebrew type: 

```bash
brew install cmake
```

- Install doxygen:

```bash
brew install doxygen
```

- Install graphViz:

```bash
brew install graphviz --with-app --with-bindings
```

- Install root6

```bash
brew install root6
```

- Install graph-tool (python). If done you can create a colored and "fancy" graph with the provided python script.

```bash
brew install graph-tool
```

- Install hdf5

    ```bash
    brew install hdf5 --with-mpi
    ```

- Install OpenMPI

    ```bash
    brew install open-mpi
    ```

- Install GSL
    
    ```bash
    brew install gsl
    ```
- Tap some repos for further sceintific packages: 

```bash
brew tap davidchall/hep
```

- Install Pythia8 for example:

```bash
brew install pythia
```

(and most of other packackes we need)

Remark: So far on brew HepMC is only available in version 2, version 3 is required for the current code, nicer wirter interfaces to root for example. So one has to install it from: http://hepmc.web.cern.ch/hepmc/

## MUSIC support

MUSIC is a (3+1)D viscous hydrodynamical code developed at McGill university.
(Official website: http://www.physics.mcgill.ca/MUSIC)
MUSIC can be integrated into the JETSCAPE framework. To download the lastest
version of MUSIC, one can run the shell script under the external_packages folder,

```bash
    ./get_music.sh
```

This shell script will clone the latest version of MUSIC to external_packages folder.
It also setup the enviroment variables for MUSIC to run. Specifically, MUSIC
needs the folder path for the EOS tables. Please make sure the enviroment
variable HYDROPROGRAMPATH to be set to the path for MUSIC code package.

When compiling MUSIC with JETSCAPE, please turn on the MUSIC support option
when generating the cmake configuration file,

```bash
    mkdir build
    cd build
    cmake -Dmusic=ON ..
    make
```

To run JETSCAPE with MUSIC, one needs to use MPI commands,

```bash
    mpirun -np 1 ./MUSICTest
```
