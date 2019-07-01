# JETSCAPE 2.0
The Jet Energy-loss Tomography with a Statistically and Computationally Advanced Program Envelope (JETSCAPE) collaboration

New with release 2.0

• TRENTO initial condition from (2+1)D to (3+1)D.

• Liquifier module to create causal jet source terms for jet energy deposition and medium evolution.

• Handling photons generated in hard scattering and in energy-loss modules. 

• Handling heavy-quarks in energy-loss modules.

• Hybrid hadronization of partons. 


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
If you are using GNU C compiler, please make sure the version of the compiler is 6 or above. 

To run JETSCAPE with MUSIC,

```bash
    ./MUSIC_evo
```

MUSIC implements parallelization using OpenMP. You can use multiple threads by setting the following enviroment variable, 

```bash
    export NUM_OMP_THREADS=4
```

This enables the MUSIC code to use 4 threads for computation.


## Using 2-stage hydro

The executable TwoStageHydro implements hydrodynamic response to jet energy-momentum deposition. 
In ./TwoStageHydro, MUSIC is run first to provide a space-time medium profile for jet shower propagation. In default, 5 hard-scattering events will be run on this profile. A new hydrodynamic background will be genereated if more than 5 hard-scattering events are computed. For each hard-scattering event, a set of source terms based on the deposited energy is generated.
MUSIC with the same initial condition is rerun with this set of source terms, to produce jets with hydrodynamic response.
To compile ./TwoStageHydro, in the build directory,

```bash
    cmake -Dmusic=ON -DiSS=ON ..
    make
```    



## SMASH hadronic afterburner

SMASH is a hadronic transport developed at Frankfurt Institute for Advanced
 Studies by the group of Prof. H. Petersen (Elfner). In JetScape SMASH can
serve as an afterburner, useful to compute soft observables.

### Installing SMASH

SMASH is published on github at https://github.com/smash-transport/smash.
See SMASH Readme for libraries required by SMASH and how to install them.

#### Prerequisites

SMASH is known to compile and work with one of these compilers (which have the
required C++11 features):
- gcc >= 4.8
- clang >= 3.2

It requires the following tools & libraries:
- cmake >= 2.8.11
- the GNU Scientific Library >= 1.15
- the Eigen3 library for linear algebra (see http://eigen.tuxfamily.org)
- boost filesystem >= 1.49
- Pythia = 8.235

See more details in SMASH README.

#### Installing Eigen

```bash
export EIGEN_DOWNLOAD_DIR=$HOME/Software
export EIGEN_INSTALL_DIR=$HOME/eigen_install

mkdir ${EIGEN_DOWNLOAD_DIR}
cd ${EIGEN_DOWNLOAD_DIR}
wget http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
tar -xf 3.2.10.tar.gz

mkdir ${EIGEN_INSTALL_DIR}
cd ${EIGEN_INSTALL_DIR}
cmake ${EIGEN_DOWNLOAD_DIR} -DCMAKE_INSTALL_PREFIX=${EIGEN_INSTALL_DIR}
make install

export EIGEN3_ROOT=${EIGEN_INSTALL_DIR}/include/eigen3/
```

Add the last export to your .bashrc file.


#### Using a custom GSL build

This is only necessary if GSL is not installed already or something
does not work with the installed version.

Download and unpack GSL:

```bash
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
    tar -zxvf gsl-latest.tar.gz
```

This creates a folder named `gsl-[version_number]` called `$GSL` here.

```bash
    cd $GSL
    ./configure --prefix $GSL
    make -jN
    make install
```

Add this export to your .bashrc file:
```bash
    export GSL_ROOT_DIR=/opt/apps/intel18/gsl/2.3
```

#### Using boost library

Assuming that boost is already installed in $HOME:

```bash
  export BOOST_ROOT=$HOME/boost_1_64_0/
```

#### Compiling SMASH library

```bash
  export JETSCAPE_DIR=${HOME}/JETSCAPE-COMP
  export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code
  export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8235
  export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8235

  cd ${JETSCAPE_DIR}/external_packages
  ./get_smash.sh
```

### Compiling JetScape with SMASH

The usage of SMASH in JetScape as an afterburner requires hydro,
sampler and SMASH itself. Therefore, to use it in JetScape,

```bash
    mkdir ${JETSCAPE_DIR}/build
    cd ${JETSCAPE_DIR}/build
    cmake -Dmusic=ON -DiSS=ON -Dfreestream=ON -Dsmash=ON ..
```

To run JetScape test with SMASH:

```bash
    cd build
    ./SMASHTest
```

Currently the iSS sampler performs resonance decays after sampling.
For reasonable physics with SMASH these decays should be switched off.
