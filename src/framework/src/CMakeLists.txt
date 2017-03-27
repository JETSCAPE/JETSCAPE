cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScape)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib/")

set (SOURCES
JetScape.cc JetScapeTask.cc
tinyxml2.cc FluidDynamics.cc HardProcess.cc JetScapeLogger.cc JetScapeXML.cc
JetEnergyLoss.cc  JetScapeModuleBase.cc JetEnergyLossManager.cc
JetScapeSignalManager.cc JetScapeWriterAscii.cc
../jet/JetClass.cpp ../jet/fjcore.cc ../jet/PartonShower.cc
PartonShowerGenerator.cc) #ElossModulesTest.cc 

if(${ZLIB_FOUND})
set (SOURCES ${SOURCES}  gzstream.cc JetScapeWriterAsciiGZ.cc)
endif()

if(${HEPMC_FOUND})
set (SOURCES ${SOURCES} JetScapeWriterHepMC.cc)
endif()

add_library(JetScape SHARED ${SOURCES})

if(${ZLIB_FOUND})
target_link_libraries(JetScape ${ZLIB_LIBRARIES})
endif()

if(${HEPMC_FOUND})
target_link_libraries(JetScape -L${HEPMC_LIB_DIR} -lHepMC)
endif()

if (${ROOT_FOUND})
target_link_libraries(JetScape ${ROOT_LIBRARIES})
endif()

#if (${PYTHIA8_FOUND})
#target_link_libraries(JetScape ${PYTHIA8_LIBRARIES})
#endif()

#target_link_libraries(JetScape JetScapeHydro GTL)
target_link_libraries(JetScape JetScapeHydro GTL)
#install(TARGETS JetScape DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

#add_library(JetScape STATIC ${SOURCES})
#Set the location for library installation -- i.e., /usr/lib in this case
#install(TARGETS xxxxx DESTINATION xxxx/lib)

# -------------------------------------------------------------------------
