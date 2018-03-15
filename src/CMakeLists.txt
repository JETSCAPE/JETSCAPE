cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScape)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "./lib")

# framework
FILE(GLOB frameworksources framework/*.cc)
set (SOURCES  ${SOURCES}  ${frameworksources} )

## Take care of options
if(NOT ${ZLIB_FOUND})
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/framework/JetScapeWriterAsciiGZ.cc)
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
if(NOT ${freestream})
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/preequilibrium/freestream-milne_jetscape.cc)
endif()

#hydro
FILE(GLOB hydrosources hydro/*.cc)
set (SOURCES  ${SOURCES}  ${hydrosources} )
## Take care of options
if(NOT ${music})
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/hydro/music_jetscape.cc)
endif()

#jet
FILE(GLOB jetsources jet/*.cc)
set (SOURCES  ${SOURCES}  ${jetsources} )

#hadronization
FILE(GLOB hadrosources hadronization/*.cc)
set (SOURCES  ${SOURCES}  ${hadrosources} )

## Take care of options
if(NOT ${iSS})
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/hadronization/iSS_jetscape.cc)
endif()

#reader
FILE(GLOB readersources reader/*.cc)
set (SOURCES  ${SOURCES}  ${readersources} )


#######################################

add_library(JetScape SHARED ${SOURCES})
target_link_libraries(JetScape JetScapeThird GTL ${PYTHIA8_LIBRARIES} libtrento ${Boost_LIBRARIES}  )

if (${ROOT_FOUND})
  target_link_libraries(JetScape ${ROOT_LIBRARIES})
endif()

if (${HDF5_FOUND})
  set_target_properties (JetScape PROPERTIES COMPILE_FLAGS "-DUSE_HDF5")
  target_link_libraries (JetScape hydroFromFile ${_hdf5_libs})
endif()

if (${HEPMC_FOUND})
  target_link_libraries (JetScape -L${HEPMC_LIB_DIR} -lHepMC )
endif()

if ( freestream )
  target_link_libraries (JetScape freestream_lib )
endif()

if ( music )
  target_link_libraries (JetScape music_lib )
endif()

if ( iSS )
  target_link_libraries (JetScape iSS_lib )
endif()
# -------------------------------------------------------------------------
