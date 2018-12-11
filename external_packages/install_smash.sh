#!/bin/bash
  module load intel/18.0.2 boost gsl cmake
  export GSL_ROOT_DIR=$(gsl-config --prefix)
  echo $GSL_ROOT_DIR

  export PYTHIAINSTALLDIR=${HOME}
  export JETSCAPE_DIR=${HOME}/JETSCAPE-COMP
  export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code

  cd ${SMASH_DIR}
  mkdir build
  cd build
  export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8230
  export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8230

  cmake .. -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
  export number_of_cores=`nproc --all`
  make -j${number_of_cores}
