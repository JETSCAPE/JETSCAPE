cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScapeTest)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib/")

set (SOURCES
     Gubser_hydro_jetscape.cc
     brick_jetscape.cc
     hydro_file_jetscape.cc
     ElossModulesTestMatter.cc
     ElossModulesTestMartini.cc
     AdSCFT.cc
     PGun.cc
     HadronizationModuleTest.cc)

if (music)
    set (SOURCES ${SOURCES} music_jetscape.cc)
endif (music)

if (freestream)
    set (SOURCES ${SOURCES} freestream-milne_jetscape.cc)
endif (freestream)

if(${PYTHIA8_FOUND})
	set (SOURCES ${SOURCES} PythiaGun.cpp)
endif()

add_library(JetScapeTest SHARED ${SOURCES})

if (${HDF5_FOUND})
    set_target_properties (JetScapeTest PROPERTIES COMPILE_FLAGS "-DUSE_HDF5")
    target_link_libraries (JetScapeTest ${_hdf5_libs})
endif()

if(${MPI_FOUND} AND ${GSL_FOUND})
    target_link_libraries (JetScapeTest ${MPI_LIBRARIES} ${GSL_LIBRARIES})
endif()

if (${PYTHIA8_FOUND})
   target_link_libraries(JetScapeTest ${PYTHIA8_LIBRARIES})
endif()

target_link_libraries(JetScapeTest JetScape JetScapeHydro hydroFromFile)

if (music)
    target_link_libraries(JetScapeTest JetScape JetScapeHydro hydroFromFile music_lib)
endif (music)

if (freestream)
    target_link_libraries(JetScapeTest JetScape JetScapeHydro hydroFromFile freestream_lib)
endif (freestream)
