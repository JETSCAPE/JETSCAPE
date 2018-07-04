cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScape)


find_library(HDF5_CXX_LIB hdf5_cpp)
if("${HDF5_CXX_LIB}" STREQUAL "HDF5_CXX_LIB-NOTFOUND")
  message(STATUS "HDF5 C++ library not found")
else()
  message(STATUS "Found HDF5 C++ library: ${HDF5_CXX_LIB}")
  # Having found the C++ lib, also locate the header and C lib.
  find_path(HDF5_INCLUDE_DIRS H5Cpp.h)
  include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
  find_library(HDF5_C_LIB hdf5)
  # Also link to zlib for HDF5 compression.
  set(HDF5_LIBRARIES ${HDF5_CXX_LIB} ${HDF5_C_LIB} z)
  # Set preprocessor definition so HDF5 output is compiled.
  add_definitions(-DTRENTO_HDF5)
endif()

if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
message(STATUS "-------------" "${HDF5_LIBRARIES}")

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "./lib")

# framework
FILE(GLOB frameworksources framework/*.cc)
set (SOURCES  ${SOURCES}  ${frameworksources} )

## Take care of options
## Note that cmake NOT needs the "" to deal with empty variables
# if(NOT "${ZLIB_FOUND}")
#   list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/framework/JetScapeWriterAsciiGZ.cc)
# endif()

if(NOT "${HEPMC_FOUND}")
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/framework/JetScapeWriterHepMC.cc)
endif()

#initialstate
FILE(GLOB initialstatesources initialstate/*.cc)
set (SOURCES  ${SOURCES}  ${initialstatesources} )

#trento needs data tables
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data_table)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/initialstate/data_table/
  DESTINATION ${CMAKE_BINARY_DIR}/data_table
  )

#pre-equilibrium
FILE(GLOB preeqsources preequilibrium/*.cc)
set (SOURCES  ${SOURCES}  ${preeqsources} )
## Take care of options
if(NOT "${freestream}")
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/preequilibrium/FreestreamMilneWrapper.cc)
endif()

#hydro
FILE(GLOB hydrosources hydro/*.cc)
set (SOURCES  ${SOURCES}  ${hydrosources} )
## Take care of options
if(NOT "${music}")
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/hydro/MusicWrapper.cc)
endif()

#jet
FILE(GLOB jetsources jet/*.cc)
set (SOURCES  ${SOURCES}  ${jetsources} )

#hadronization
FILE(GLOB hadrosources hadronization/*.cc)
set (SOURCES  ${SOURCES}  ${hadrosources} )

## Take care of options
if(NOT "${iSS}")
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/hadronization/iSpectraSamplerWrapper.cc)
endif()

#reader
FILE(GLOB readersources reader/*.cc)
set (SOURCES  ${SOURCES}  ${readersources} )


#######################################

add_library(JetScape SHARED ${SOURCES})
target_link_libraries(JetScape JetScapeThird GTL ${PYTHIA8_LIBRARIES} libtrento libLBTD ${HDF5_LIBRARIES} ${GSL_LIBRARIES} ${GSLCALAS_LIBRARIES} ${Boost_LIBRARIES} ${_hdf5_libs} -pthread -lpthread   )

if (${ROOT_FOUND})
  target_link_libraries(JetScape ${ROOT_LIBRARIES})
endif()

if (${HDF5_FOUND})
  target_link_libraries (JetScape hydroFromFile ${_hdf5_libs})
endif()

if ( ${HEPMC_FOUND} )
  target_link_libraries (JetScape -L${HEPMC_LIB_DIR} -lHepMC )
endif()

if ( music )
  target_link_libraries (JetScape music_lib )
endif()

if ( iSS )
  target_link_libraries (JetScape iSS_lib )
endif()
# -------------------------------------------------------------------------
