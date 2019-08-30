module load cmake/3.7.1
module load intel/18.0.2
module load gsl
module load boost
module load hdf5
module load impi
module load eigen

#export EIGEN_DOWNLOAD_DIR=$HOME/Software
export EIGEN_INSTALL_DIR=$TACC_EIGEN_DIR
export EIGEN3_ROOT=$TACC_EIGEN_INC
export GSL=$(gsl-config --prefix)
export GSL_HOME=$(gsl-config --prefix)
export GSL_ROOT_DIR=$(gsl-config --prefix)

export JETSCAPE_DIR=/work/05780/everett/stampede2/JETSCAPE-COMP-RUN-EBE

export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code
export number_of_cores=`nproc --all`

export PYTHIAINSTALLDIR=$HOME

export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8235
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8235

export CC=icc
export CXX=icpc
export OpenMP_CXX=icpc
