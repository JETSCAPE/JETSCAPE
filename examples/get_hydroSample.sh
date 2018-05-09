#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
 *
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

# download VISHNU hydro examples
curlcmd=wget
command -v ${curlcmd} > /dev/null || curlcmd="curl -LO"
command -v ${curlcmd} > /dev/null || { echo "Please install curl or wget" ; exit 1; }
$curlcmd "https://bitbucket.org/sscao/hydro_profiles/downloads/PbPb2760-Avg-00-05.tar.gz"

tar xvzf PbPb2760-Avg-00-05.tar.gz
mkdir -p test_hydro_files
mv PbPb2760-Avg-00-05 test_hydro_files/event-0
rm PbPb2760-Avg-00-05.tar.gz
