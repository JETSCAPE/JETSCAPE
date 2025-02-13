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

# default LHAPATH environment variable
export LHAPATH='/home/jetscape-user/.local/share/LHAPDF'

# create the LHAPATH if it doesn't exist
if [ ! -d "$LHAPATH" ]; then
    mkdir -p $LHAPATH
fi

# replace lhapdf_set with the desired set or pass as an argument
lhapdf_set = "JAM20-SIDIS_PDF_proton_nlo"

if [ $# -eq 1 ]; then
    lhapdf_set=$1
fi

wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/$lhapdf_set.tar.gz -O- | tar xz -C $LHAPATH
