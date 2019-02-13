############################################################
#
# This script runs a unit test which runs TRENTO_FS_HYDRO and SAMPLE_AFTERBURNER
#
# DO NOT FORGET TO CHANGE THE EMAIL ADDRESS IN THE MPI JOB SCRIPT
# Written by Matthew Heffernan, based off readme by Weiyao Ke in sims branch commit a909ee0
# and slack messages from Weiyao and Derek Everett
###########################################################################################

module load intel/18.0.2 cmake/3.7.1 gsl boost hdf5 eigen impi # Load libraries on stampede2

WORKDIR=$WORK #$(pwd) # see sims_install.sh for discussion on how/why to change working directory
cd $WORKDIR
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

# Create test directory and populate with the necessary files
cd $JETSCAPE
/bin/mkdir -p workflow_test_unit && cp -r build/. workflow_test_unit

cd workflow_test_unit

# A few additional files are needed from sims_scripts
cp ../workflow_test/generate_module_input_files.py .
cp ../workflow_test/Run_WorkflowTest.py .
cp ../workflow_test/workflow_test .
cp ../workflow_test/workflowtest_oversampled_afterburner.py .
cp -r ../workflow_test/smash_input/ .

# Ensuring the proper input data/configurations are in place
cp -r ../external_packages/iS3D/PDG/ .
cp -r ../external_packages/iS3D/tables/ .
cp -r ../external_packages/iS3D/deltaf_coefficients/ .
cp -r ../external_packages/music/EOS/ .

# Generate the input files and submit a job to stampede2 
python generate_module_input_files.py
sbatch workflow_test
