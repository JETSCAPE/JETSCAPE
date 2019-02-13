###########################################################################################
# Compiles the sims branch in the current working directory on stampede2.tacc.utexas.edu
# This is for preparation to run a unit test if the user has changed anything
#
#
# Written by Matthew Heffernan, based off readme by Weiyao Ke in sims branch commit a909ee0
###########################################################################################

WORKDIR=$WORK #$(pwd) # to install the sims branch in the current directory, uncomment here and change mpi_job_skylake submission script accordingly
cd $WORKDIR

#Setup environment for Pythia
export PATH=$WORKDIR/software/bin:$PATH
export LD_LIBRARY_PATH=$WORKDIR/software/lib:\$LD_LIBRARY_PATH
export PYTHIA8DATA=$WORKDIR/software/share/Pythia8/xmldoc/
export PYTHIA8=$WORKDIR/software

##Compilation on stampede2 TACC

# Prepare environment
JETSCAPE=$WORKDIR/JETSCAPE-COMP

module load cmake/3.7.1 # Include these again since setting environment variables seems to reset these
module load intel/18.0.2
module load gsl
module load boost
module load hdf5
module load impi
module load eigen

export EIGEN_DOWNLOAD_DIR=$HOME/Software
export EIGEN_INSTALL_DIR=$TACC_EIGEN_DIR
export EIGEN3_ROOT=$TACC_EIGEN_INC
export GSL=$(gsl-config --prefix)
export GSL_HOME=$(gsl-config --prefix)
export GSL_ROOT_DIR=$(gsl-config --prefix)
export JETSCAPE_DIR=$WORKDIR/JETSCAPE-COMP
export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code
export number_of_cores=1 #$(nproc) #`nproc --all`
export PYTHIAINSTALLDIR=$WORKDIR
export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8235
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8235

export CC=icc
export CXX=icpc
export OpenMP_CXX=icpc

# Compile and install JETSCAPE
cd $JETSCAPE
pwd
/bin/mkdir -p build && cd build
pwd
cmake -Dmusic=on -DiS3D=on -Dsmash=on -Dfreestream=on ../

make
make install
