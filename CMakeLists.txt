# ##############################################################################
# Initial Cmake Setup ###
# ##############################################################################

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project (JetScape CXX C)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message("source directory:" ${CMAKE_SOURCE_DIR})
  message("binary directory:" ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# Default to the release build.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# for mac compliance
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0054 NEW)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmakemodules")

# ##############################################################################
# Parse options       ###
# ##############################################################################

# Unit Tests. Turn off with 'cmake -Dunittests=OFF'.
option(unittests "Build all unittests." ON)

# freestream. Turn on with 'cmake -DUSE_FREESTREAM=ON'. Note that some warnings
# are generated. Could be turned off by adding the following to CFLAGS in
# external_packages/freestream-milne/Makefile -Wno-unknown-pragmas
# -Wno-writable-strings -Wno-return-type
# -Wc++11-compat-deprecated-writable-strings
option(USE_FREESTREAM "Build tests for freestream-milne" OFF)
if(USE_FREESTREAM)
  message("Includes for freestream ...")
  include_directories(./external_packages/freestream-milne/src)
endif(USE_FREESTREAM)

# CLVisc. Trun on with 'cmake -DUSE_CLVISC=ON'.
find_package(OpenCL)

option(USE_CLVISC "Build tests for CLVisc" OFF)
if(OPENCL_FOUND)
  message("Found OpenCL, try clvisc with cmake -DUSE_CLVISC=ON")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
  if(USE_CLVISC)
    message("Includes for clvisc ...")
    include_directories(./external_packages/clvisc_wrapper/include/)
    include_directories(./external_packages/clvisc_wrapper/)
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external_packages/PyVisc")
      message(
        FATAL_ERROR
          "Error: CLVisc source has not been downloaded in external_packages by ./get_clvisc.sh"
      )
    endif()
  endif(USE_CLVISC)
endif(OPENCL_FOUND)

# IP-Glasma. Turn on with 'cmake -DUSE_IPGLASMA=ON'.
option(USE_IPGLASMA "Build with IP-Glasma" OFF)
if(USE_IPGLASMA)
  message("Includes for IP-Glasma ...")
  include_directories(./external_packages/ipglasma
                      ./external_packages/ipglasma/src)
endif(USE_IPGLASMA)

# MUSIC. Turn on with 'cmake -DUSE_MUSIC=ON'.
option(USE_MUSIC "Build tests for MUSIC" OFF)
if(USE_MUSIC)
  message("Includes for MUSIC ...")
  include_directories(./external_packages/music ./external_packages/music/src)
endif(USE_MUSIC)

# Soft Particlization. Turn on with 'cmake -DUSE_MUSIC=ON -DUSE_ISS=ON'.
option(USE_ISS "Build tests for iSS" OFF)
if(USE_ISS)
  message("Includes for iSS ...")
  include_directories(./external_packages/iSS ./external_packages/iSS/src)
endif(USE_ISS)

# SMASH afterburner. Turn on with 'cmake -DUSE_MUSIC=ON -DUSE_ISS=ON
# -DSMASH=ON'.
option(USE_SMASH "Build tests for SMASH" OFF)
if(USE_SMASH)
  message("SMASH includes and library ...")
  find_package(SMASH)
  if(${SMASH_FOUND})
    include_directories(${SMASH_INCLUDE_DIR})
  endif(${SMASH_FOUND})
endif(USE_SMASH)

# Compile with OpenMP support
if(APPLE)
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
  set(OpenMP_CXX_LIB_NAMES "libomp")
  set(OpenMP_libomp_LIBRARY "/usr/local/opt/libomp/lib")
endif()
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# ##############################################################################
# Compiler & Linker Flags ###
# ##############################################################################
message("Compiler and Linker flags ...")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -pipe -Wall")
# can turn off some warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder -Wno-unused-variable ")
# Then set the build type specific options. These will be automatically
# appended.
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# can turn on debugging information set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# system dependence. note that APPLE also defines UNIX, hence the elseif to
# differentiate
if(APPLE)
  message(STATUS "Apple : " ${CMAKE_HOST_SYSTEM})
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-but-set-variable -Wno-parentheses -fext-numeric-literals"
    )
  endif()

  # can turn off some warnings set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}
  # -Wno-unused-private-field")
elseif(UNIX)
  message(STATUS "Linux : " ${CMAKE_HOST_SYSTEM})
  # can turn off some warnings
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-but-set-variable -Wno-parentheses -fext-numeric-literals"
    )
  endif()

  # Additional options set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs
  # -ftest-coverage")
endif()

message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
# message(STATUS "LD_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})

# ##############################################################################
# add the library path and inclusion path of trento to jetscape Placement of
# this package is important. Needs to come after CXX flags (cause it needs to
# pickup fPIC) and before install options (cause it overwrites the default
# location)
add_subdirectory(external_packages/trento/)
include_directories("${CMAKE_SOURCE_DIR}/external_packages/trento/src/")
# ##############################################################################

add_subdirectory(external_packages/googletest/)
include_directories(
  "${CMAKE_SOURCE_DIR}/external_packages/googletest/googletest/include/")

# ##############################################################################
# Installing Header and Library Files   ###
# ##############################################################################
# Run with, e.g., cmake -DCMAKE_INSTALL_PREFIX=~/tmp .. default directory is the
# build directory Note that trento also automatically populates a bin directory
# and we cannot disable this behavior Also, iSS and mpihydro install binaries in
# CMAKE_HOME_DIR. Sigh.
#
# Install header files default install prefix: build directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
      "${PROJECT_BINARY_DIR}"
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
    "${CMAKE_SOURCE_DIR}/src/afterburner/"
    "${CMAKE_SOURCE_DIR}/src/jet/"
    "${CMAKE_SOURCE_DIR}/src/reader/"
    # external packages
    "${CMAKE_SOURCE_DIR}/external_packages/"
    "${CMAKE_SOURCE_DIR}/external_packages/ipglasma/src/"
    "${CMAKE_SOURCE_DIR}/external_packages/iSS/src/"
    "${CMAKE_SOURCE_DIR}/external_packages/hydro_from_external_file/src/"
    "${CMAKE_SOURCE_DIR}/external_packages/music/src/"
    "${CMAKE_SOURCE_DIR}/external_packages/trento/src/"
    "${CMAKE_SOURCE_DIR}/external_packages/clvisc_wrapper/"
    "${CMAKE_SOURCE_DIR}/external_packages/smash/src/"
    "${CMAKE_SOURCE_DIR}/external_packages/googletest/googletest/src/"
  DESTINATION "include" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "*.h*" # select header files
  # Necessary to exclude directories to prevent a whole (mostly empty) hierarchy
  PATTERN "gtl" EXCLUDE
  PATTERN "iSS" EXCLUDE
  PATTERN "hydro_from_external_file" EXCLUDE
  PATTERN "ipglasma" EXCLUDE
  PATTERN "music" EXCLUDE
  PATTERN "clvisc" EXCLUDE
  PATTERN "trento" EXCLUDE
  PATTERN "tests" EXCLUDE
  PATTERN "data_table" EXCLUDE
  PATTERN "LBT-tables" EXCLUDE
  PATTERN "Martini" EXCLUDE
  PATTERN "googletest" EXCLUDE)

# We have includes of the form #include "GTL/dijkstra.h" which needs to be
# handled separately
install(
  # source directory
  DIRECTORY "${CMAKE_SOURCE_DIR}/external_packages/gtl/include/"
  DESTINATION "include" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "GTL/*.h*" # select header files
  # Necessary to exclude directories to prevent a whole (mostly empty) hierarchy
  # PATTERN "gtl/src" EXCLUDE
)

install(
  # source directory
  DIRECTORY "${CMAKE_SOURCE_DIR}/external_packages/clvisc_wrapper/include/"
  DESTINATION "include" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "*.h*" # select header files
  # Necessary to exclude directories to prevent a whole (mostly empty) hierarchy
)

# Install lib files
install(
  # our libraries
  DIRECTORY "${PROJECT_BINARY_DIR}/src/lib/"
            "${PROJECT_BINARY_DIR}/lib/"
            # external packages
            "${PROJECT_BINARY_DIR}/external_packages/gtl/lib/"
            "${PROJECT_BINARY_DIR}/external_packages/ipglasma/src/"
            "${PROJECT_BINARY_DIR}/external_packages/iSS/src/"
            "${PROJECT_BINARY_DIR}/external_packages/music/src/"
            "${PROJECT_BINARY_DIR}/external_packages/trento/src/"
            "${PROJECT_BINARY_DIR}/external_packages/clvisc_wrapper/"
  DESTINATION "lib" # target directory
  FILES_MATCHING # install only matched files
  PATTERN "lib*.*" # selects .so, .a, .dylib, ...
  PATTERN "CMakeFiles" EXCLUDE)

# ##############################################################################
# Required packages   ###
# ##############################################################################

# Find and use Boost.
message("Looking for Boost ...")
find_package(Boost 1.50 REQUIRED COMPONENTS filesystem program_options system)
message(STATUS "Boost Include dirs : " ${Boost_INCLUDE_DIRS})
# boost needs special treatment. For reasons.
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

message("Looking for ZLIB ...")
find_package(ZLIB REQUIRED)
message(STATUS "ZLib found")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_GZIP")

message("Looking for Pythia8 ...")
find_package(Pythia8 REQUIRED)
include_directories(${PYTHIA8_INCLUDE_DIR})
message(STATUS "Pythia8 Include dir : " ${PYTHIA8_INCLUDE_DIR})

message("Looking for HDF5 ...")
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HDF5")
set(_hdf5_libs ${HDF5_CXX_LIB} ${HDF5_C_LIB} ${ZLIB} ${HDF5_HL_LIBRARIES})
message(
  STATUS "HDF5 include dir and libs : ${HDF5_INCLUDE_DIRS} ${_hdf5_libs} ")

# ##############################################################################
# Optional packages   ###
# ##############################################################################
unset(USE_HEPMC)
message("Looking for HepMC ...")
find_package(HEPMC)
if(${HEPMC_FOUND})
  include_directories(${HEPMC_INCLUDE_DIR})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HEPMC")
  message(STATUS "HepMC Include dir : " ${HEPMC_INCLUDE_DIR})
endif()

option(USE_ROOT "Build using ROOT Libraries and Output" OFF)
if(USE_ROOT)
  message("Looking for ROOT ...")
  # find_package(ROOT)
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
  find_package(ROOT)
  if(${ROOT_FOUND})
    # include_directories(${ROOT_INCLUDES})
    include_directories(${ROOT_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_ROOT")
    message(STATUS "ROOT Include dir : " ${ROOT_INCLUDE_DIRS})
    message(STATUS "ROOT Libraries : " ${ROOT_LIBRARIES})
  endif()
endif(USE_ROOT)

message("Looking for MPI ...")
find_package(MPI)
if(${MPI_FOUND})
  message(STATUS "MPI Include dir : " ${MPI_INCLUDE_PATH})
  include_directories(${MPI_INCLUDE_PATH})
endif(${MPI_FOUND})

message("Looking for GSL ...")
find_package(GSL)
if(${GSL_FOUND})
  message(STATUS "GSL_INCLUDE_DIR : " ${GSL_INCLUDE_DIR}
                 ${GSLCBLAS_INCLUDE_DIRS})
  include_directories(${GSL_INCLUDE_DIR} ${GSLCBLAS_INCLUDE_DIRS})
endif(${GSL_FOUND})

# ##############################################################################
# Framework Includes  ###
# ##############################################################################
message("Include Directories ...")
include_directories(./src/)
include_directories(./src/framework)
include_directories(./src/initialstate)
include_directories(./src/preequilibrium)
include_directories(./src/liquefier)
include_directories(./src/hydro)
include_directories(./src/hadronization)
include_directories(./src/afterburner)
include_directories(./src/jet)
include_directories(./src/reader)
include_directories(./external_packages/)
include_directories(./external_packages/gtl/include)
include_directories(./external_packages/hydro_from_external_file)
include_directories(./external_packages/Random123)
include_directories(./examples/)
if(APPLE)
  # set OpenMP_INCLUDE_DIRS to the path of the OpenMP headers for example:
  # /usr/local/opt/libomp/include
  include_directories($ENV{OpenMP_INCLUDE_DIRS})
endif()
if("${ROOT_FOUND}")
  include_directories(./src/root)
  # add_subdirectory(./src/root)
endif()

# include_directories will be expanded as packages are found

# ##############################################################################
# Some additional settings for subdirectories ###
# ##############################################################################
add_subdirectory(./external_packages)
add_subdirectory(./external_packages/gtl)
add_subdirectory(./src)

if(unittests)
  add_subdirectory(./examples/unittests/)
endif(unittests)

if(USE_SMASH)
  target_compile_definitions(JetScape PRIVATE USE_SMASH)
endif(USE_SMASH)

if(USE_FREESTREAM)
  if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external_packages/freestream-milne")
    message(
      FATAL_ERROR
        "Error: freestream-milne source has not been downloaded in external_packages by ./get_freestream-milne.sh"
    )
  endif()
  message("Building freestream-milne ...")
  if(${GSL_FOUND})
    add_subdirectory(./external_packages/freestream-milne)
  else()
    message(FATAL_ERROR "GSL are necessary for Freestreaming")
  endif()
  target_compile_definitions(JetScape PRIVATE USE_FREESTREAM)
endif(USE_FREESTREAM)

if (USE_IPGLASMA)
  if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external_packages/ipglasma")
      message(FATAL_ERROR "Error: IP-Glasma source has not been downloaded in external_packages by ./get_ipglasma.sh")
  endif()
  if (${GSL_FOUND})
    message("Building IP-Glasma ...")
    add_subdirectory(./external_packages/ipglasma)
  else()
    message (FATAL_ERROR "GSL is necessary for IP-Glasma" )
  endif()
  target_compile_definitions(JetScape PRIVATE USE_IPGLASMA)
endif()

if(USE_MUSIC)
  if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external_packages/music")
    message(
      FATAL_ERROR
        "Error: MUSIC source has not been downloaded in external_packages by ./get_music.sh"
    )
  endif()
  if(${GSL_FOUND})
    message("Building MUSIC ...")
    add_subdirectory(./external_packages/music)
  else()
    message(FATAL_ERROR "GSL are necessary for MUSIC")
  endif()
  target_compile_definitions(JetScape PRIVATE USE_MUSIC)
endif()

if(USE_ISS)
  if(EXISTS "${CMAKE_SOURCE_DIR}/external_packages/iSS")
    add_subdirectory(${CMAKE_SOURCE_DIR}/external_packages/iSS)
  else()
    message(
      FATAL_ERROR
        "Error: iSS source has not been downloaded in external_packages by ./get_iSS.sh"
    )
  endif()
  target_compile_definitions(JetScape PRIVATE iSpectraSampler)
endif(USE_ISS)

if(OPENCL_FOUND AND USE_CLVISC)
  add_subdirectory(${CMAKE_SOURCE_DIR}/external_packages/clvisc_wrapper)
endif(OPENCL_FOUND AND USE_CLVISC)

if(OPENMP_FOUND)
  target_compile_definitions(JetScape PRIVATE USE_OPENMP)
endif()

# ##############################################################################
# Binary location     ###
# ##############################################################################
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# ##############################################################################
# Standalone Reader lib  ###
# ##############################################################################
file(GLOB LIBREADERSOURCES src/reader/*.cc)
set(LIBREADERSOURCES ${LIBREADERSOURCES})
set(LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/JetScapeParticles.cc)
set(LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/StringTokenizer.cc)
set(LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/JetClass.cc)
set(LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/JetScapeLogger.cc)
set(LIBREADERSOURCES ${LIBREADERSOURCES} src/framework/PartonShower.cc)

add_library(JetScapeReader SHARED ${LIBREADERSOURCES})
set_target_properties(JetScapeReader PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                                ${PROJECT_BINARY_DIR}/lib)
target_link_libraries(JetScapeReader JetScapeThird GTL ${PYTHIA8_LIBRARIES})

# ##############################################################################
# Executables         ###
# ##############################################################################

# Run Jetscape
add_executable(runJetscape ./examples/runJetscape.cc)
target_link_libraries(runJetscape JetScape)

# Read Jetscape output
add_executable(readerTest ./examples/readerTest.cc)
target_link_libraries(readerTest JetScape)

add_executable(FinalStateHadrons ./examples/FinalStateHadrons.cc)
target_link_libraries(FinalStateHadrons JetScape)

add_executable(FinalStatePartons ./examples/FinalStatePartons.cc)
target_link_libraries(FinalStatePartons JetScape)

# executables with additional dependencies
if ( USE_IPGLASMA )
    target_link_libraries(runJetscape ${GSL_LIBRARIES})
endif (USE_IPGLASMA )
if ( USE_MUSIC AND USE_ISS )
    target_link_libraries(runJetscape ${GSL_LIBRARIES})
endif (USE_MUSIC AND USE_ISS )

if(USE_MUSIC
   AND USE_ISS
   AND USE_SMASH)
  include_directories(${SMASH_INCLUDE_DIR})
  # "-lubsan -lasan")
  target_link_libraries(runJetscape ${SMASH_LIBRARIES} ${MPI_LIBRARIES}
                        ${GSL_LIBRARIES})
endif(
  USE_MUSIC
  AND USE_ISS
  AND USE_SMASH)

# if (${ROOT_FOUND}) add_executable(PythiaBrickTestRoot
# ./examples/custom_examples/PythiaBrickTestRoot.cc)
# target_link_libraries(PythiaBrickTestRoot JetScape) endif(${ROOT_FOUND})

# -----------------------------------------------------------
#### copy essential files for MUSIC to build/
file(COPY ./external_packages/trento/nucleusConfigs DESTINATION ${CMAKE_BINARY_DIR}/)
if (USE_IPGLASMA)
  if(EXISTS "${CMAKE_SOURCE_DIR}/external_packages/ipglasma/")
    file(COPY ${CMAKE_SOURCE_DIR}/examples/test_ipglasma_files/ipglasma.input
         DESTINATION ${CMAKE_BINARY_DIR}/)
    file(
      COPY ${CMAKE_SOURCE_DIR}/external_packages/ipglasma/qs2Adj_vs_Tp_vs_Y_200.in
      DESTINATION ${CMAKE_BINARY_DIR}/)
  else()
    message(FATAL_ERROR "Error: Cannot find files for IP-Glasma.")
  endif()
endif()
if(USE_FREESTREAM)
  if(EXISTS "${CMAKE_SOURCE_DIR}/examples/test_freestream_files/")
    file(COPY ./examples/test_freestream_files/freestream_input
         DESTINATION ${CMAKE_BINARY_DIR})
  else()
    message(FATAL_ERROR "Error: Cannot find files for freestream-milne.")
  endif()
endif(USE_FREESTREAM)
if(USE_MUSIC)
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

# copy essential files for iSS to build/
if(USE_ISS)
  if(EXISTS "${CMAKE_SOURCE_DIR}/external_packages/iSS/iSS_tables")
    file(COPY ${CMAKE_SOURCE_DIR}/external_packages/iSS/iSS_tables
         DESTINATION ${CMAKE_BINARY_DIR}/)
    file(COPY ${CMAKE_SOURCE_DIR}/external_packages/iSS/iSS_parameters.dat
         DESTINATION ${CMAKE_BINARY_DIR}/)
  else()
    message(
      FATAL_ERROR
        "Error: iSS tables have not been downloaded in external_packages by ./get_iSS.sh"
    )
  endif()
endif(USE_ISS)

# copy essential files for LBT to build/
if(EXISTS "${CMAKE_SOURCE_DIR}/external_packages/LBT-tables")
  execute_process(
    COMMAND
      ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_SOURCE_DIR}/external_packages/LBT-tables
      ${CMAKE_BINARY_DIR}/LBT-tables)
else()
  message(
    "Warning: LBT-tables have not been downloaded in external_packages by ./get_lbtTab.sh; LBT-brickTest will NOT run properly."
  )
endif()
