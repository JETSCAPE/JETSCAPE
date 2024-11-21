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

# download the 3+1D OSU freestreaming code
#git clone --depth=1 https://github.com/derekeverett/freestream-milne.git freestream-milne
#git clone --depth=1 https://github.com/derekeverett/freestream-milne.git -b time_step_history freestream-milne
#git clone --depth=1 https://github.com/chunshen1987/freestream-milne -b time_step_history freestream-milne

# using a commit from the freestream-milne repository that is compatible with JETSCAPE 3.6.6
folderName="freestream-milne"
#commitHash="94722958595cb712fdb00cc59375ad7c9030faed"
commitHash="0d980593396a8b2f6a9b3933780215df00022d58"

#git clone https://github.com/chunshen1987/freestream-milne -b time_step_history $folderName
git clone https://github.com/chunshen1987/freestream-milne -b fix_output_orientation $folderName
cd $folderName
git checkout $commitHash

#cd freestream-milne
#patch -p0 -Ni ../freestream-milne-external-params.patch
