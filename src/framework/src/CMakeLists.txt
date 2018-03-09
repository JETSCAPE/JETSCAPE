cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScape)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib/")

set (SOURCES
  JetScape.cc
  JetScapeTask.cc
  tinyxml2.cc
  FluidDynamics.cc
  PreequilibriumDynamics.cc
  HardProcess.cc
  Collision.cc
  TrentoInitial.cc
  JetScapeLogger.cc
  JetScapeXML.cc
  JetEnergyLoss.cc
  JetScapeModuleBase.cc
  JetEnergyLossManager.cc
  JetScapeSignalManager.cc
  JetScapeWriterAscii.cc
  JetScapeEvent.cc
  JetScapeTaskSupport.cc
  JetScapeParticles.cc
  ../jet/JetClass.cpp
  ../jet/fjcore.cc
  ../jet/PartonShower.cc
  InitialState.cc
  PartonShowerGenerator.cc
  PartonPrinter.cc
  HadronizationManager.cc
  Hadronization.cc
  SoftParticlization.cc
) 

if(${ZLIB_FOUND})
    set (SOURCES ${SOURCES}  gzstream.cc JetScapeWriterAsciiGZ.cc)
endif()

if (${HEPMC_FOUND})
    set (SOURCES ${SOURCES} JetScapeWriterHepMC.cc)
endif(${HEPMC_FOUND})

add_library(JetScape SHARED ${SOURCES})

if(${ZLIB_FOUND})
    target_link_libraries(JetScape ${ZLIB_LIBRARIES})
endif()

if (${HEPMC_FOUND})
    set_target_properties(JetScape PROPERTIES COMPILE_FLAGS "-DUSE_HEPMC")
    target_link_libraries(JetScape -L${HEPMC_LIB_DIR} -lHepMC)
endif()

if (${ROOT_FOUND})
    target_link_libraries(JetScape ${ROOT_LIBRARIES})
endif()

add_subdirectory(test)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data_table)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data_table/
            DESTINATION ${CMAKE_BINARY_DIR}/data_table)

target_link_libraries(JetScape libtrento ${Boost_LIBRARIES} ${HDF5_LIBRARIES})

#target_link_libraries(JetScape JetScapeHydro GTL)
target_link_libraries(JetScape JetScapeHydro GTL ${PYTHIA8_LIBRARIES} )
#install(TARGETS JetScape DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

#add_library(JetScape STATIC ${SOURCES})
#Set the location for library installation -- i.e., /usr/lib in this case
#install(TARGETS xxxxx DESTINATION xxxx/lib)

# -------------------------------------------------------------------------
