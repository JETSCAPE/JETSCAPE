cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (testJetScape)

cmake_policy(SET CMP0042 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/Modules")

message("Looking for ZLIB ...")
find_package(ZLIB)
if(${ZLIB_FOUND})
   message(STATUS "ZLib found")
endif()

message("Looking for HepMC ...")
find_package(HEPMC)

message("Looking for ROOT ...")
find_package(ROOT)
    
message("Looking for Pythia8 ...")
find_package(Pythia8)

message("Looking for HDF5 ...")
find_package(HDF5)
if(${HDF5_FOUND})
    include_directories(${HDF5_INCLUDE_DIRS})
    set(_hdf5_libs ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif(${HDF5_FOUND})

message("Looking for MPI ...")
find_package(MPI)
if(${MPI_FOUND})
    include_directories(${MPI_INCLUDE_PATH})
endif(${MPI_FOUND})

message("Looking for GSL ...")
find_package(GSL)
if(${GSL_FOUND})
    include_directories(${GSL_INCLUDE_DIR})
endif(${GSL_FOUND})


include_directories(${GSL_INCLUDE_DIR})

message ("Include Directories ...")

include_directories(./src ./test ./hydro ./jet ./gtl/include ./reader) #for test case ...
include_directories(./3rdparty/hydro_from_external_file)
include_directories(./3rdparty/music ./3rdparty/music/src)

if (${HEPMC_FOUND})
      include_directories(${HEPMC_INCLUDE_DIR})
      Message(STATUS "HepMC Include dir : " ${HEPMC_INCLUDE_DIR})
endif()

if (${PYTHIA8_FOUND})
      include_directories(${PYTHIA8_INCLUDE_DIR})
     Message(STATUS "Pythia8 Include dir : " ${PYTHIA8_INCLUDE_DIR})
endif()

if (${ROOT_FOUND})
    include_directories(${ROOT_INCLUDES})
    Message(STATUS "ROOT Include dir : " ${ROOT_INCLUDES})
    #Message(STATUS ${ROOT_LIBRARIES})
    
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    #set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -m64")    
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./bin/")

message("Compiler and Linker flags ...")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -fPIC -pipe -Wall -g -Wno-unused-variable -Wno-unused-private-field -Wno-reorder")

if(LINUX) 
    message( STATUS "Linux : " ${CMAKE_HOST_SYSTEM})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if (APPLE)
     message( STATUS "Apple : " ${CMAKE_HOST_SYSTEM})
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-macro-redefined -std=c++11 -stdlib=libc++")
     if(${MPI_FOUND})
          set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
          set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
     endif(${MPI_FOUND})
     message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})	
     message(STATUS "LD_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})
endif()   

# w/o not compiling library automatically
add_subdirectory(./gtl)
add_subdirectory(./src)
add_subdirectory(./hydro) #for test case fluid_dynamics ...
add_subdirectory(./jet)  #for test case jet classes …
add_subdirectory(./test)  #for test case ...
add_subdirectory(./reader)  #for test case ...
add_subdirectory(./3rdparty/hydro_from_external_file)  #for test case ...
add_subdirectory(./3rdparty/music)  #for test case ...

link_directories(./lib)

add_executable(testJetScape ./src/testJetScape.cc)
add_executable(brickTest ./src/brickTest.cc)
add_executable(readerTest ./reader/readerTest.cc)
add_executable(hydroFileTest ./src/hydroFileTest.cc)
add_executable(MUSICTest ./src/MUSICTest.cc)

target_link_libraries(testJetScape JetScape JetScapeHydro JetScapeTest hydroFromFile MUSIC) #for test case ...
target_link_libraries(brickTest JetScape JetScapeHydro  JetScapeTest) #for test case ...
target_link_libraries(readerTest JetScape JetScapeHydro  JetScapeReader) #for reader test case ...
target_link_libraries(hydroFileTest JetScape JetScapeHydro JetScapeTest hydroFromFile) #for test case ...
target_link_libraries(MUSICTest JetScape JetScapeHydro JetScapeTest MUSIC) #for test case ...

if (${HDF5_FOUND})
    set_target_properties (hydroFileTest PROPERTIES COMPILE_FLAGS "-DUSE_HDF5")
    target_link_libraries (hydroFileTest ${_hdf5_libs})
endif()

if (${MPI_FOUND} AND ${GSL_FOUND})
    target_link_libraries (MUSICTest ${MPI_LIBRARIES} ${GSL_LIBRARIES})
endif()

#Check how to set env variables ...
if (${ROOT_FOUND})
   target_link_libraries(testJetScape ${ROOT_LIBRARIES})
   target_link_libraries(brickTest ${ROOT_LIBRARIES})
   target_link_libraries(hydroFileTest ${ROOT_LIBRARIES})
endif()

if (${PYTHIA8_FOUND})
   target_link_libraries(testJetScape ${PYTHIA8_LIBRARIES})
   target_link_libraries(brickTest ${PYTHIA8_LIBRARIES})
   target_link_libraries(hydroFileTest ${PYTHIA8_LIBRARIES})
   message( ${PYTHIA8_LIBRARIES} )
endif()

#add_executable(testJS ./src/testJetScape.cc)
#target_link_libraries(testJS JetScape)

# -----------------------------------------------------------
    
#set(LIBS ${LIBS} ${ROOT_LIBRARIES})
#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp")
