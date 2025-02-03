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

# using a commit from the freestream-milne repository that is compatible with the current JETSCAPE version
folderName="freestream-milne"
commitHash="e0a21feb48a922b4b4541ab0e4d745c65594bb5f"

git clone https://github.com/chunshen1987/freestream-milne -b time_step_history $folderName
cd $folderName
git checkout $commitHash
