#!/bin/sh

if [ -z $LD_LIBRARY_PATH ]; then
export LD_LIBRARY_PATH
fi

if [ -z $DYLD_LIBRARY_PATH ]; then 
export DYLD_LIBRARY_PATH
fi

export BASEDIR=${HOME}

export JetScape=${PWD}/lib
export LD_LIBRARY_PATH=${JetScape}:${LD_LIBRARY_PATH}
#only for Mac needed
export DYLD_LIBRARY_PATH=${JetScape}:${DYLD_LIBRARY_PATH}

# PYTHIA8 directory
export PYTHIA8DIR=${BASEDIR}/pythia8223
export PYTHIA8_ROOT_DIR=${BASEDIR}/pythia8223
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
