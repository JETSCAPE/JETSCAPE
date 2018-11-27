#export env variables
export EIGEN_DOWNLOAD_DIR=$HOME/Software
export EIGEN_INSTALL_DIR=$HOME/eigen_install
export EIGEN3_ROOT=${EIGEN_INSTALL_DIR}/include/eigen3/
export GSL=/home1/05780/everett/gsl-2.5
export GSL_HOME=/home1/05780/everett/gsl-2.5
export GSL_ROOT_DIR=/opt/apps/intel18/gsl/2.3
export JETSCAPE_DIR=${WORK}/JETSCAPE-COMP-SIMS
export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code
export number_of_cores=`nproc --all`
export PYTHIAINSTALLDIR=/home1/05780/everett/install_pythia_for_smash
export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8230
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8230

#load necessary modules
module load cmake/3.7.1
module load intel/18.0.0
module load boost
module load hdf5
module load impi

#set compiler
export CC=icc
export CXX=icpc
export OpenMP_CXX=icpc

#for gsl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home1/05780/everett/gsl-2.5/lib
