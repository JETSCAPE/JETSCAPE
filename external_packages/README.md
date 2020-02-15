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
    cmake -DUSE_MUSIC=ON ..
    make
```

To run JETSCAPE with MUSIC, one needs to use MPI commands,

```bash
    mpirun -np 1 ./MUSICTest
```

## Running JETSCAPE with CLVisc
In order to run clvisc in jetscape, one has to download it in external\_packages/, using 

```
sh get_clvisc.sh
```

Then compile and run the framework with a XML configuration file which turns clvisc on.
```
cd build/
cmake .. -DUSE_CLVISC=on
make
./runJetscape ../config/jetscape_clvisc.xml
```
If the cmake fails because OpenCL is not installed, please check it.
OpenCL is delivered in MacOS by default. 
If you use linux machine with Nvidia GPUs, you will need to install CUDA,
which will provide OpenCL support.
If you use linux machine with AMD GPUs, Intel GPUs or any CPUs,
you will need to install AMD APP SDK.

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
    cmake -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_FREESTREAM=ON -DUSE_SMASH=ON ..
```

To run JetScape test with SMASH:

```bash
    cd build
    ./SMASHTest
```

Currently the iSS sampler performs resonance decays after sampling.
For reasonable physics with SMASH these decays should be switched off.
