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

#
# 1) Download the SMASH code
#    Currently it only works if one has an access to smash repository,
#    which is not public. In future the repository has to be replaced to
#    a public one.
git clone hyihp-repos@fias.uni-frankfurt.de:scm/smash.git smash/smash_code

#
# 2) Checkout to the right branch
#    This is only necessary now, when the changes in SMASH, that allow
#    connection to JetScape, are not merged to master
#    (and are not part of the tagged version).
cd smash/smash_code
git checkout oliiny/SMASH_with_JetScape

#
# 3) Compile SMASH
#
mkdir build
cd build
cmake ..
number_of_cores=`nproc --all`
echo "Compiling SMASH using ${number_of_cores} cores."
make -j${number_of_cores} SmashShared

#
# 4) Copy experiment.cc to includes
#    Unfortunately I(oliiny) have found no other way to include templated
#    smash::Experiment<Modus> to JetScape, except dangerous
#    #include "smash/experiment.cc".
cd ..
cp src/experiment.cc src/include/smash/

#
# 5) Set the environment variable, that allows JetScape to find SMASH
#
cd ../..
export SMASH_DIR=${PWD}/smash/smash_code
echo "SMASH_DIR is set to ${SMASH_DIR}."
echo "It won't work it you used ./get_smash.sh: use \". get_smash.sh\" instead."
