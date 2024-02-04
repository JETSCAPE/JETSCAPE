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
git clone --depth=1 https://github.com/chunshen1987/freestream-milne -b time_step_history freestream-milne

#cd freestream-milne
#patch -p0 -Ni ../freestream-milne-external-params.patch
