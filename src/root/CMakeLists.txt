cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScapeRoot)

#for MacOSX
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0022 NEW)

if(${HEPMC_FOUND})

  set (SOURCES ${SOURCES} JetScapeWriterRootHepMC.cc)
  set (HEADERS ${HEADERS} JetScapeWriterRootHepMC.h)

  message( STATUS "ROOT Sources : " ${SOURCES})

  #---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
  include(${ROOT_USE_FILE})

  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib/")
  set(CMAKE_INSTALL_LIBDIR "./")

  #message("ROOT_GENERATE_DICTIONARY .. ")
  include_directories(./)
  add_definitions(${ROOT_CXX_FLAGS})

  #ROOT_GENERATE_DICTIONARY(G__JetScapeRootLibDict ${HEADERS} LINKDEF JetScapeRootLinkDef.h)

  #set (SOURCES ${SOURCES} G__JetScapeRootLibDict.cxx)

  add_library(JetScapeRoot SHARED ${SOURCES}) # G__JetScapeRootLibDict.cxx)
  #add_dependencies(JetScapeRoot G__JetScapeRootLibDict)

  target_link_libraries(JetScapeRoot JetScape ${ROOT_LIBRARIES} -L${HEPMC_LIB_DIR} -lHepMC3 -lHepMC3rootIO)

endif()
