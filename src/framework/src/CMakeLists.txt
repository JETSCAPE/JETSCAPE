cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScape)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib/")

set (SOURCES
JetScape.cc JetScapeTask.cc
tinyxml2.cc FluidDynamics.cc HardProcess.cc JetScapeLogger.cc JetScapeXML.cc
JetEnergyLoss.cc  JetScapeModuleBase.cc JetEnergyLossManager.cc
JetScapeSignalManager.cc JetScapeWriterAscii.cc) #ElossModulesTest.cc 

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

target_link_libraries(JetScape JetScapeHydro JetScapeEloss)

#install(TARGETS JetScape DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

#add_library(JetScape STATIC ${SOURCES})
#Set the location for library installation -- i.e., /usr/lib in this case
#install(TARGETS xxxxx DESTINATION xxxx/lib)

# -------------------------------------------------------------------------

#include_directories(./)

#if (HEPMC_FOUND)
#     include_directories(${HEPMC_INCLUDE_DIR}/include)
#endif()

#else (HEPMC_FOUND)
 #    include_directories(./)
#endif(HEPMC_FOUND)

#file(GLOB SOURCES "./*.cc")

#if (!${CMAKE_CXX_FLAGS})
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -fPIC -pipe -Wall -g")
#endif()

#if(LINUX) 
    #message( ">>> Linux : " ${CMAKE_HOST_SYSTEM})
        # linux stuff here
#endif()

#if (APPLE)
     #message( ">>> Apple : " ${CMAKE_HOST_SYSTEM})
     #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
     #message(">>> CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
#endif()   

#set_property(TARGET gtest APPEND_STRING PROPERTY COMPILE_FLAGS " -w")
#https://www.johnlamp.net/cmake-tutorial-6-realistically-getting-a-boost.html
#find_package(Boost 1.32 REQUIRED COMPONENTS program_options)
#include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

#if (NOT CONFIGURED_ONCE)
#    set(CMAKE_CXX_FLAGS "${warnings}"
#        CACHE STRING "Flags used by the compiler during all build types." FORCE)
#    set(CMAKE_C_FLAGS   "${warnings}"
#        CACHE STRING "Flags used by the compiler during all build types." FORCE)
#endif()

# see also : $ROOTSYS/cmake !!!!
# https://root.cern.ch/how/integrate-root-my-project-cmake
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
#find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#include(${ROOT_USE_FILE})

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#ROOT_GENERATE_DICTIONARY(G__Event Event.h LINKDEF EventLinkDef.h)

#---Create a shared library with geneated dictionary
#add_library(Event SHARED Event.cxx G__Event.cxx)
#target_link_libraries(Event ${ROOT_LIBRARIES})
