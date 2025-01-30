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

SMASH [https://smash-transport.github.io] is a hadronic transport approach
developed at Frankfurt University and GSI by the group of
Prof. H. Elfner (nee Petersen).  In JetScape SMASH can
serve as an afterburner, useful to compute soft observables.

### Installing SMASH

SMASH is published on github at https://github.com/smash-transport/smash.
See SMASH Readme for libraries required by SMASH and how to install them.

```
  export EIGEN3_ROOT=<eigen install directory>/include/eigen3/
  export GSL_ROOT_DIR=$(gsl-config --prefix)
  export BOOST_ROOT=<boost install directory>
  export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8235
  export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8235

  export JETSCAPE_DIR=${HOME}/JETSCAPE-COMP
  export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code

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

## IP-Glasma

IP-Glasma requires the use of FFTW. The package `libfftw3-dev` is already included in JETSCAPE's supported [Docker image](https://hub.docker.com/r/jetscape/base).

### Installing IP-Glasma

To download IP-Glasma, one can run the shell script under the external_packages folder.

```bash
./get_ipglasma.sh
```

This shell script will clone IP-Glasma to the external_packages folder.

### Compiling JETSCAPE with IP-Glasma

Create a build folder and cd into it.

```bash
cd ${JETSCAPE_DIR}
mkdir build
cd build
```

When compiling IP-Glasma with JETSCAPE, please turn on the IP-Glasma support option when generating the cmake configuration file.

```bash
cmake .. -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_FREESTREAM=ON -DUSE_SMASH=ON -DUSE_IPGlasma=ON

make -j4  # Builds using 4 cores; adapt as appropriate
```
If you're not compiling all the external modules, you don't have to turn them on.
