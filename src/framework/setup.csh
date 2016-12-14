#!/bin/csh

if(! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH

setenv BASEDIR ${HOME}

setenv JetScape ${PWD}/lib
setenv LD_LIBRARY_PATH ${JetScape}:${LD_LIBRARY_PATH}
#only for Mac needed
setenv DYLD_LIBRARY_PATH ${JetScape}:${DYLD_LIBRARY_PATH}

if ($?TERM == 0 || $?prompt == 0) exit 0

echo ''
echo 'Setup JetScape Library'
echo '======================'
echo ''
echo "<I>---------------Info--------------------<I>"
echo "Setting up the following environments: "
echo "JetScape: " $JetScape
echo "<I>---------------Info--------------------<I>"
echo ""
