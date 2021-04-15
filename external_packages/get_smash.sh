#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

# 1) Download the SMASH code
git clone https://github.com/smash-transport/smash.git smash/smash_code

# 2) Compile SMASH
cd smash/smash_code
#checkout the commit version that was used for production runs
git checkout 0063efcc88c11151fa4422940a8bd145a52c356d
# Need to specify a dummy name and email here for the cherry-pick to succeed...
git config user.email "jetscape@jetscape.org"
git config user.name "JETSCAPE"
# Include fixes for g++ 10
git cherry-pick 353b1b5f55c72867044927029566528b67720e9f
# Fixes for pythia version detection
git cherry-pick 5f05981d1f5f7904b184c3dfdc599718211aa6f9
# Print the log to verify that we're where we think we are.
git log -n 10 --oneline

#update March 7, 2019 Smash Bug fixed w.r.t. formation time
#if necessary for future stability, checkout a fixed commit after this date
rm -r build
mkdir build
cd build
cmake .. -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
number_of_cores=`nproc --all`
number_of_cores_to_compile=$(( ${number_of_cores} > 20 ? 20 : ${number_of_cores} ))
echo "Compiling SMASH using ${number_of_cores_to_compile} cores."
make -j${number_of_cores_to_compile} smash_shared
