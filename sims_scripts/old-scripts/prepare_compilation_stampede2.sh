#export env variables
#load necessary modules
module load cmake/3.7.1
module load intel/18.0.2
module load gsl
module load boost
module load hdf5
module load impi

export EIGEN_DOWNLOAD_DIR=$HOME/Software
export EIGEN_INSTALL_DIR=$HOME/eigen_install
export EIGEN3_ROOT=${EIGEN_INSTALL_DIR}/include/eigen3/
export GSL=$(gsl-config --prefix)
export GSL_HOME=$(gsl-config --prefix)
export GSL_ROOT_DIR=$(gsl-config --prefix)
export JETSCAPE_DIR=$WORK/fresh_start/JETSCAPE-COMP-SIMS
export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code
export number_of_cores=`nproc --all`
export PYTHIAINSTALLDIR=$HOME
export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8235
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8235

export CC=icc
export CXX=icpc
export OpenMP_CXX=icpc

