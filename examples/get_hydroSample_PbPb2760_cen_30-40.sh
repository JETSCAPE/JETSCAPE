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

if [ -d "test_hydro_files/event-0" ]; then rm -Rf test_hydro_files/event-0; fi

# download VISHNU hydro examples
curlcmd=wget
command -v ${curlcmd} > /dev/null || curlcmd="curl -LO"
command -v ${curlcmd} > /dev/null || { echo "Please install curl or wget" ; exit 1; }
$curlcmd "https://bitbucket.org/sscao/hydro_profiles/downloads/PbPb2760-Avg-30-40.tar.gz"

tar xvzf PbPb2760-Avg-30-40.tar.gz
mkdir -p test_hydro_files
mv PbPb2760-Avg-30-40 test_hydro_files/event-0
rm PbPb2760-Avg-30-40.tar.gz
