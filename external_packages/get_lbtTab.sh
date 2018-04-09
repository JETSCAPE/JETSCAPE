#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2017
#
# For the full list of contributors see AUTHORS.
# Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
# or via email to bugs.jetscape.org@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

 
# download the tables required by LBT code
curlcmd=wget
command -v ${curlcmd} > /dev/null || curlcmd="curl -LO"
command -v ${curlcmd} > /dev/null || { echo "Please install curl or wget" ; exit 1; }
$curlcmd "https://bitbucket.org/sscao/lbt-tables/downloads/LBT-tables.tar.gz"

tar xvzf LBT-tables.tar.gz

