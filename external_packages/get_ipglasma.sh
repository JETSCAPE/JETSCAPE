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

folderName="ipglasma"
commitHash="215aea40fd3e0777f01c464cc4031fb2b4344449"

# download the code package
git clone --depth=1 -b InterfaceJETSCAPE https://github.com/chunshen1987/ipglasma $folderName
cd $folderName
git checkout $commitHash