cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScape)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "./lib")

# framework
FILE(GLOB frameworksources framework/*.cc)
set (SOURCES  ${SOURCES}  ${frameworksources} )
# set (SOURCES
#   ${SOURCES}
#   framework/JetScape.cc
#   framework/JetScapeModuleBase.cc
#   framework/JetScapeReader.cc
#   framework/JetScapeSignalManager.cc
#   framework/JetScapeTask.cc
#   framework/JetScapeTaskSupport.cc
#   framework/JetScapeLogger.cc
#   framework/JetScapeParticles.cc
#   framework/JetClass.cpp
#   framework/JetScapeWriterAscii.cc
#   framework/JetScapeXML.cc
#   framework/StringTokenizer.cc
# )

## Take care of options
if(NOT ${ZLIB_FOUND})
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/framework/JetScapeWriterAsciiGZ.cc)
endif()

#initialstate
FILE(GLOB initialstatesources initialstate/*.cc)
set (SOURCES  ${SOURCES}  ${initialstatesources} )
## Take care of options
if(NOT ${iSS})
  list (REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/initialstate/iSS_jetscape.cc)
endif()

#trento needs data tables
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data_table)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/initialstate/data_table/
  DESTINATION ${CMAKE_BINARY_DIR}/data_table
  )

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

#hadro
FILE(GLOB hadrosources hadro/*.cc)
set (SOURCES  ${SOURCES}  ${hadrosources} )

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
  target_link_libraries (JetScape hydroFromFile ${_hdf5_libs})
endif()

if (${HEPMC_FOUND})
  target_link_libraries (JetScape -L${HEPMC_LIB_DIR} -lHepMC )
endif()

if ( music )
  target_link_libraries (JetScape music_lib )
endif()

if ( iSS )
  target_link_libraries (JetScape iSS_lib )
endif()
# -------------------------------------------------------------------------
