###########################
### Initial Cmake Setup ###
###########################

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (JetScape)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message("source directory:" ${CMAKE_SOURCE_DIR})
  message("binary directory:" ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# for mac compliance
cmake_policy(SET CMP0042 NEW)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmakemodules")

###########################
### Parse options       ###
###########################

# Unit Tests. Turn on with 'cmake -Dtest=ON'.
option(test "Build all tests." OFF)

# freestream. Turn on with 'cmake -Dfreestream=ON'.
# Note that some warnings are generated. Could be turned off by adding the following to CFLAGS in
# external_packages/freestream-milne/Makefile
# -Wno-unknown-pragmas -Wno-writable-strings -Wno-return-type -Wc++11-compat-deprecated-writable-strings
option(freestream "Build tests for freestream-milne" OFF)
if (freestream)
  message("Includes for freestream ...")
  include_directories(./external_packages/freestream-milne/src)
endif (freestream)

# MUSIC. Turn on with 'cmake -Dmusic=ON'.
option(music "Build tests for MUSIC" OFF)
if (music)
    message("Includes for music ...")
    include_directories(./external_packages/music ./external_packages/music/src)
    find_package(OpenMP REQUIRED)
endif (music)

# Soft Particlization. Turn on with 'cmake -Dmusic=ON -DiSS=ON'.
option(iSS "Build tests for iSS" OFF)
if (iSS)
  message("Includes for iSS ...")
  include_directories(./external_packages/iSS ./external_packages/iSS/src)
endif (iSS)

###############################
### Compiler & Linker Flags ###
###############################
message("Compiler and Linker flags ...")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -fPIC -pipe -Wall -std=c++11")
## can turn off some warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder -Wno-unused-variable ")
## can turn on debugging information
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

### system dependence.
### note that APPLE also defines UNIX, hence the elseif to differentiate
if (APPLE)
    message( STATUS "Apple : " ${CMAKE_HOST_SYSTEM})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") 
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-but-set-variable -Wno-parentheses")
    endif()

    ## can turn off some warnings
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-private-field")
elseif(UNIX) 
    message( STATUS "Linux : " ${CMAKE_HOST_SYSTEM})
    ## can turn off some warnings
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-but-set-variable -Wno-parentheses")
    endif()

    ## Additional options
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()   

message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})	
#message(STATUS "LD_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})

###########################
# add the library path and inclusion path of trento to jetscape
# Placement of this package is important.
# Needs to come after CXX flags (cause it needs to pickup fPIC)
# and before install options (cause it overwrites the default location)
add_subdirectory(external_packages/trento/)
include_directories("${CMAKE_SOURCE_DIR}/external_packages/trento/src/")
###########################

#############################################
### Installing Header and Library Files   ###
#############################################
## Run with, e.g., cmake -DCMAKE_INSTALL_PREFIX=~/tmp .. 
## default directory is the build directory
## Note that trento also automatically populates a bin directory
## and we cannot disable this behavior
## Also, iSS and mpihydro install binaries in CMAKE_HOME_DIR. Sigh.
###
# Install header files
# default install prefix: build directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}"
    CACHE PATH "Install path prefix, prepended onto install directories."
    FORCE)
endif()
message(STATUS "Now: install prefix is ${CMAKE_INSTALL_PREFIX}")

install(
  # source directory
  DIRECTORY
  # our source
  "${CMAKE_SOURCE_DIR}/src/framework/"
  "${CMAKE_SOURCE_DIR}/src/hadronization/"
  "${CMAKE_SOURCE_DIR}/src/initialstate/"
  "${CMAKE_SOURCE_DIR}/src/hydro/"
  "${CMAKE_SOURCE_DIR}/src/jet/"
  "${CMAKE_SOURCE_DIR}/src/reader/"
  # external packages
  "${CMAKE_SOURCE_DIR}/external_packages/"
  "${CMAKE_SOURCE_DIR}/external_packages/iSS/src/"
  "${CMAKE_SOURCE_DIR}/external_packages/hydro_from_external_file/src/"
  "${CMAKE_SOURCE_DIR}/external_packages/music/src/"
  "${CMAKE_SOURCE_DIR}/external_packages/trento/src/"
  DESTINATION "include" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "*.h*" # select header files
  ## Necessary to exclude directories to prevent a whole (mostly empty) hierarchy
  PATTERN "gtl" EXCLUDE
  PATTERN "iSS" EXCLUDE
  PATTERN "hydro_from_external_file" EXCLUDE
  PATTERN "music" EXCLUDE
  PATTERN "trento" EXCLUDE
  PATTERN "tests" EXCLUDE
  PATTERN "data_table" EXCLUDE
  PATTERN "LBT-tables" EXCLUDE
  PATTERN "Martini" EXCLUDE
  PATTERN "googletest" EXCLUDE
  )

## We have includes of the form #include "GTL/dijkstra.h"
## which needs to be handled separately
install(
  # source directory
  DIRECTORY
  "${CMAKE_SOURCE_DIR}/external_packages/gtl/include/"
  DESTINATION "include" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "GTL/*.h*" # select header files
  ## Necessary to exclude directories to prevent a whole (mostly empty) hierarchy
  # PATTERN "gtl/src" EXCLUDE
  )

# Install lib files
install(
  # our libraries
  DIRECTORY
  "${PROJECT_BINARY_DIR}/src/lib/"
  "${PROJECT_BINARY_DIR}/lib/"
  # external packages
  "${PROJECT_BINARY_DIR}/external_packages/gtl/lib/"
  "${PROJECT_BINARY_DIR}/external_packages/iSS/src/"
  "${PROJECT_BINARY_DIR}/external_packages/music/src/"
  "${PROJECT_BINARY_DIR}/external_packages/trento/src/"
  DESTINATION "lib" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "lib*.*" # selects .so, .a, .dylib, ...
  PATTERN "CMakeFiles" EXCLUDE
  )

###########################
### Required packages   ###
###########################


# Find and use Boost.
message("Looking for Boost ...")
find_package(Boost 1.50 REQUIRED COMPONENTS filesystem program_options system)
Message(STATUS "Boost Include dirs : " ${Boost_INCLUDE_DIRS})
## boost needs special treatment. For reasons.
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

message("Looking for ZLIB ...")
find_package(ZLIB)
if(${ZLIB_FOUND})
  message(STATUS "ZLib found")
  set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DUSE_GZIP" )
endif()

###########################
### Optional packages   ###
###########################
message("Looking for HepMC ...")
find_package(HEPMC)
if (${HEPMC_FOUND})
  include_directories(${HEPMC_INCLUDE_DIR})
  set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DUSE_HEPMC" )
  Message(STATUS "HepMC Include dir : " ${HEPMC_INCLUDE_DIR})
endif()

message("Looking for ROOT ...")
find_package(ROOT)
if (${ROOT_FOUND})
  include_directories(${ROOT_INCLUDES})
  Message(STATUS "ROOT Include dir : " ${ROOT_INCLUDES})
  Message(STATUS "ROOT Libraries : " ${ROOT_LIBRARIES})
endif()

message("Looking for Pythia8 ...")
find_package(Pythia8)
if (${PYTHIA8_FOUND})
  include_directories(${PYTHIA8_INCLUDE_DIR})
  Message(STATUS "Pythia8 Include dir : " ${PYTHIA8_INCLUDE_DIR})
endif()

message("Looking for HDF5 ...")
find_package(HDF5)
if(${HDF5_FOUND})
  include_directories(${HDF5_INCLUDE_DIRS})
  message(STATUS "HDF5 Include dir : " ${HDF5_INCLUDE_DIRS})
  set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DUSE_HDF5" )
  set(_hdf5_libs ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif(${HDF5_FOUND})

message("Looking for MPI ...")
find_package(MPI)
if(${MPI_FOUND})
  message(STATUS "MPI Include dir : " ${MPI_INCLUDE_PATH})
  include_directories(${MPI_INCLUDE_PATH})
endif(${MPI_FOUND})

message("Looking for OpenMP ...")
find_package(OpenMP)
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )


message("Looking for GSL ...")
find_package(GSL)
if(${GSL_FOUND})
  message(STATUS "GSL_INCLUDE_DIR : " ${GSL_INCLUDE_DIR})
  include_directories(${GSL_INCLUDE_DIR})
endif(${GSL_FOUND})


###########################
### Framework Includes  ###
###########################
message ("Include Directories ...")
include_directories(./src/ )
include_directories(./src/framework )
include_directories(./src/initialstate )
include_directories(./src/preequilibrium )
include_directories(./src/hydro )
include_directories(./src/hadronization )
include_directories(./src/jet )
include_directories(./src/reader )
include_directories(./external_packages/)
include_directories(./external_packages/gtl/include )
include_directories(./external_packages/hydro_from_external_file )
include_directories(./examples/ )


### include_directories will be expanded as packages are found

###################################################
### Some additional settings for subdirectories ###
###################################################
add_subdirectory(./external_packages)
add_subdirectory(./external_packages/gtl)
add_subdirectory(./src)

if (freestream)
  if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external_packages/freestream-milne")
    message(FATAL_ERROR "Error: freestream-milne source has not been downloaded in external_packages by ./get_freestream-milne.sh")
  endif()
  message("Building freestream-milne ...")
  if (${GSL_FOUND})
    add_subdirectory(./external_packages/freestream-milne/src)
  else()
      message (FATAL_ERROR "GSL are necessary for Freestreaming" )
  endif()
endif (freestream)

if (music)
  if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external_packages/music")
    message(FATAL_ERROR "Error: MUSIC source has not been downloaded in external_packages by ./get_music.sh")
  endif()  
  if (${OpenMP_CXX_FOUND} AND ${GSL_FOUND})
    message("Building MUSIC ...")
    add_subdirectory(./external_packages/music)
  else()
    message (FATAL_ERROR "OpenMP and GSL are necessary for MUSIC" )
  endif()
endif()

if (iSS)
  if(EXISTS "${CMAKE_SOURCE_DIR}/external_packages/iSS")
    add_subdirectory( ${CMAKE_SOURCE_DIR}/external_packages/iSS )
  else()
    message(FATAL_ERROR "Error: iSS source has not been downloaded in external_packages by ./get_iSS.sh")
  endif()      
endif (iSS)



###########################
### Binary location     ###
###########################
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})


##############################
### Standalone Reader lib  ###
##############################
FILE(GLOB LIBREADERSOURCES src/reader/*.cc)
set (LIBREADERSOURCES ${LIBREADERSOURCES} )
set (LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/JetScapeParticles.cc )
set (LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/StringTokenizer.cc )
set (LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/JetClass.cc )
set (LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/JetScapeLogger.cc )
set (LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/PartonShower.cc )

add_library(JetScapeReader SHARED ${LIBREADERSOURCES})
set_target_properties(JetScapeReader PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
target_link_libraries(JetScapeReader JetScapeThird GTL ${PYTHIA8_LIBRARIES})

###########################
### Executables         ###
###########################

### Run Jetscape
add_executable(brickTest ./examples/brickTest.cc)
target_link_libraries(brickTest JetScape )

add_executable(PythiaBrickTest ./examples/PythiaBrickTest.cc)
target_link_libraries(PythiaBrickTest JetScape )

add_executable(LBT_brickTest ./examples/LBT_brickTest.cc)
target_link_libraries(LBT_brickTest JetScape )

if(${HDF5_FOUND})
  add_executable(hydroFileTest ./examples/hydroFileTest.cc)
  target_link_libraries(hydroFileTest JetScape ${_hdf5_libs})
  
  add_executable(hydroJetTest  ./examples/hydroJetTest.cc)
  target_link_libraries(hydroJetTest JetScape ${_hdf5_libs})

  add_executable(hydroJetTestPGun  ./examples/hydroJetTestPGun.cc)
  target_link_libraries(hydroJetTestPGun JetScape ${_hdf5_libs})
endif()

add_executable(SimpleValidate ./examples/SimpleValidate.cc)
target_link_libraries(SimpleValidate JetScape )

### Read Jetscape output
add_executable(readerTest ./examples/readerTest.cc)
target_link_libraries(readerTest JetScape )

add_executable(FinalStateHadrons ./examples/FinalStateHadrons.cc)
target_link_libraries(FinalStateHadrons JetScape )

add_executable(FinalStatePartons ./examples/FinalStatePartons.cc)
target_link_libraries(FinalStatePartons JetScape )

# executables with additional dependencies
if (freestream AND music)
  add_executable(freestream-milneTest ./examples/freestream-milneTest.cc)
  target_link_libraries(freestream-milneTest JetScape ${GSL_LIBRARIES})
endif (freestream AND music)
if ( music AND iSS )
  add_executable(MUSICTest ./examples/MUSICTest.cc)
  target_link_libraries(MUSICTest JetScape ${OpenMP_CXX_LIBRARIES} ${GSL_LIBRARIES})
endif (music AND iSS )

# if (${ROOT_FOUND})
#   add_executable(pwg2_reader ./examples/pwg2_reader.cxx)
#   target_link_libraries(pwg2_reader JetScape )
# endif(${ROOT_FOUND})

# if (${ROOT_FOUND})
#   target_link_libraries(brickTest ${ROOT_LIBRARIES})
#   target_link_libraries(LBT_brickTest ${ROOT_LIBRARIES})
#   target_link_libraries(PythiaBrickTest ${ROOT_LIBRARIES})
#   target_link_libraries(writerTest ${ROOT_LIBRARIES})
# endif()


# -----------------------------------------------------------
#### copy jetscape_init.xml to build/, such that one can test xml reading in build/ directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples/jetscape_init.xml DESTINATION ${CMAKE_BINARY_DIR}/)

#### copy essential files for MUSIC to build/
if (freestream)
  if(EXISTS "${CMAKE_SOURCE_DIR}/examples/test_freestream_files/")
    file(COPY ./examples/test_freestream_files/freestream_input
      DESTINATION ${CMAKE_BINARY_DIR})
  else()
    message(FATAL_ERROR "Error: Cannot find files for freestream-milne.")
  endif()
endif (freestream)
if (music)
  if(EXISTS "${CMAKE_SOURCE_DIR}/external_packages/music/EOS/")
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/EOS)
    file(COPY ${CMAKE_SOURCE_DIR}/external_packages/music/EOS/
      DESTINATION ${CMAKE_BINARY_DIR}/EOS)
    file(COPY ${CMAKE_SOURCE_DIR}/examples/test_music_files/music_input
      DESTINATION ${CMAKE_BINARY_DIR})
  else()
    message(FATAL_ERROR "Error: Cannot find files for MUSIC.")
  endif()
endif()

#### copy essential files for iSS to build/
if (iSS)
  if(EXISTS "${CMAKE_SOURCE_DIR}/external_packages/iSS/iSS_tables")
    file(COPY ${CMAKE_SOURCE_DIR}/external_packages/iSS/iSS_tables DESTINATION ${CMAKE_BINARY_DIR}/)
    file(COPY ${CMAKE_SOURCE_DIR}/external_packages/iSS/iSS_parameters.dat DESTINATION ${CMAKE_BINARY_DIR}/)
  else()
    message(FATAL_ERROR "Error: iSS tables have not been downloaded in external_packages by ./get_iSS.sh")
  endif()
endif (iSS)

#### copy essential files for LBT to build/
if(EXISTS "${CMAKE_SOURCE_DIR}/external_packages/LBT-tables")
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/external_packages/LBT-tables ${CMAKE_BINARY_DIR}/LBT-tables)
else()
  message("Warning: LBT-tables have not been downloaded in external_packages by ./get_lbtTab.sh; LBT-brickTest will NOT run properly.")
endif()

