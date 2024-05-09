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

# download the code package
# git clone --depth=1 https://github.com/chunshen1987/iSS -b JETSCAPE iSS

# using a commit from the iSS repository that is compatible with JETSCAPE 3.6.1
folderName="iSS"
commitHash="9d1847c8b6e02be6f018cf5b1003fdf6085d328f"

git clone https://github.com/chunshen1987/iSS -b JETSCAPE $folderName
cd $folderName
git checkout $commitHash
