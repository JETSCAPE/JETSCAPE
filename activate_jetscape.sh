#!/bin/sh

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

if [ -z $LD_LIBRARY_PATH ]; then
export LD_LIBRARY_PATH
fi

if [ -z $DYLD_LIBRARY_PATH ]; then 
export DYLD_LIBRARY_PATH
fi

export BASEDIR=${HOME}
export PYTHIAINSTALLDIR=${HOME}

export JetScape=${PWD}/lib
export LD_LIBRARY_PATH=${JetScape}:${LD_LIBRARY_PATH}
#only for Mac needed
export DYLD_LIBRARY_PATH=${JetScape}:${DYLD_LIBRARY_PATH}

# PYTHIA8 directory
export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8235
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8235
#export PYTHIA8_INCLUDE_DIR=`${PYTHIA8DIR}/bin/pythia8-config --includedir`/Pythia8
#export PYTHIA8_LIBRARIES=`${PYTHIA8DIR}/bin/pythia8-config --libdir`
export LD_LIBRARY_PATH=${PYTHIA8DIR}/lib:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${PYTHIA8DIR}/lib:${DYLD_LIBRARY_PATH}

#ROOT setup
export ROOTSYS="/opt/local/libexec/root5"
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

if [ -z ${TERM} -o -z ${SHELL} ]; then exit 0
fi

echo ''
echo 'Setup JetScape Library'
echo '======================'
echo ''
echo "<I>---------------Info--------------------<I>"
echo "Setting up the following environments: "
echo "JetScape: " ${JetScape}
echo "Pythia8: "${PYTHIA8DIR}/lib    
echo "<I>---------------Info--------------------<I>"
echo ""
