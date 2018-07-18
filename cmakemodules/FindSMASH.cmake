# Try to find SMASH instalation
#
# Once done this will define:
#
# SMASH_LIBRARIES      smash libraries
# SMASH_INCLUDE_DIR    directory of smash includes
# SMASH_FOUND          true is smash found
#
# The environment variable SMASH_DIR must be set properly to succeed, e. g.:
# export SMASH_DIR=~/Work/SMASH/smash

message(STATUS "Looking for SMASH ...")

set(SMASH_INCLUDE_DIR
   $ENV{SMASH_DIR}/3rdparty/Cuba-4.2
   $ENV{SMASH_DIR}/3rdparty/einhard
   $ENV{SMASH_DIR}/3rdparty/yaml-cpp-0.6.2/include
   $ENV{SMASH_DIR}/3rdparty/pythia8230/include
   $ENV{SMASH_DIR}/src/include
)
message(STATUS "SMASH includes found in ${SMASH_INCLUDE_DIR}")

#find_library(SMASH_LIBRARY NAMES smash_lib PATHS $ENV{SMASH_DIR}/build/src)
#find_library(EINHARD_LIBRARY NAMES einhard PATHS $ENV{SMASH_DIR}/build/3rdparty/einhard)
#find_library(CPPYAML_LIBRARY NAMES yaml-cpp PATHS $ENV{SMASH_DIR}/build/3rdparty/yaml-cpp-0.6.2)
#find_library(SMASH_PYTHIA_LIBRARY NAMES pythia PATHS $ENV{SMASH_DIR}/build/3rdparty/pythia8230)
#find_library(INTEGRATION_LIBRARY NAMES cuhre PATHS $ENV{SMASH_DIR}/build/3rdparty/Cuba-4.2/src/cuhre)
#set(SMASH_LIBRARIES ${EINHARD_LIBRARY} ${CPPYAML_LIBRARY} ${SMASH_PYTHIA_LIBRARY} ${SMASH_LIBRARY} ${INTEGRATION_LIBRARY})

find_library(SMASH_LIBRARIES NAMES SmashShared PATHS $ENV{SMASH_DIR}/build/src)

message(STATUS "SMASH libraries: ${SMASH_LIBRARIES}")

set(SMASH_FOUND FALSE)
if (SMASH_INCLUDE_DIR AND SMASH_LIBRARIES)
#    SMASH_LIBRARY     AND
#    EINHARD_LIBRARY   AND
#    CPPYAML_LIBRARY   AND
#    SMASH_PYTHIA_LIBRARY AND
#    INTEGRATION_LIBRARY)
   set(SMASH_FOUND TRUE)
endif ()

message(STATUS "SMASH found: ${SMASH_FOUND}")
