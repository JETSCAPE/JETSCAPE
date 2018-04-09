#!/bin/csh

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

if(! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH

setenv BASEDIR ${HOME}

setenv JetScape ${PWD}/lib
setenv LD_LIBRARY_PATH ${JetScape}:${LD_LIBRARY_PATH}
#only for Mac needed
setenv DYLD_LIBRARY_PATH ${JetScape}:${DYLD_LIBRARY_PATH}

# PYTHIA8 directory
setenv PYTHIA8DIR ${BASEDIR}/pythia8223
setenv LD_LIBRARY_PATH ${PYTHIA8DIR}/lib:${LD_LIBRARY_PATH}
setenv DYLD_LIBRARY_PATH ${PYTHIA8DIR}/lib:${DYLD_LIBRARY_PATH}
    
if ($?TERM == 0 || $?prompt == 0) exit 0

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
