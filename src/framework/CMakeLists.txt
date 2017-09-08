cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (testJetScape)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message("source directory:" ${CMAKE_SOURCE_DIR})
    message("binary directory:" ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

cmake_policy(SET CMP0042 NEW)

# Tell cmake where to find modules
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

message ("Include Directories ...")

include_directories(./src ./test ./hydro ./jet ./gtl/include ./reader) #for test case ...
include_directories(./3rdparty/hydro_from_external_file)

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

if (music)
    message("Includes for music ...")
    find_package(MPI)
    if(${MPI_FOUND})
        include_directories(${MPI_INCLUDE_PATH})
    endif(${MPI_FOUND})
    message("Looking for GSL ...")
    find_package(GSL)
    if(${GSL_FOUND})
        include_directories(${GSL_INCLUDE_DIR})
    endif(${GSL_FOUND})
    include_directories(./3rdparty/music ./3rdparty/music/src)
endif (music)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./bin/")

message("Compiler and Linker flags ...")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -fPIC -pipe -Wall -g -Wno-reorder -Wno-unused-variable ")

if (APPLE)
     message( STATUS "Apple : " ${CMAKE_HOST_SYSTEM})

     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-private-field -Wno-macro-redefined -std=c++11 -stdlib=libc++")
     if(${MPI_FOUND})
          set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
          set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
     endif(${MPI_FOUND})

     message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})	
     message(STATUS "LD_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})
elseif(UNIX) 
    message( STATUS "Linux : " ${CMAKE_HOST_SYSTEM})
    #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -fprofile-arcs -ftest-coverage -fPIC -pipe -Wall -g -Wno-reorder -Wno-unused-variable") 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-sign-compare -Wno-unused-but-set-variable -Wno-parentheses -fprofile-arcs -ftest-coverage")
    #set(CMAKE_BUILD_TYPE Debug)
endif()   


# have binary in the build directory, easier to access
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# Find and use Boost.
find_package(Boost 1.50 REQUIRED COMPONENTS filesystem program_options system)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})


# add the library path and inclusion path of trento to jetscape
add_subdirectory(3rdparty/trento/)
include_directories("${CMAKE_SOURCE_DIR}/3rdparty/trento/src/")

# Options. Turn on with 'cmake -Dtest=ON'.
option(test "Build all tests." ON) # Makes boolean 'test' available.

# Options. Turn on with 'cmake -Dmusic=ON'.
option(music "Build tests for MUSIC" OFF) # Makes boolean 'test' available.

# w/o not compiling library automatically
add_subdirectory(./gtl)
add_subdirectory(./src)
add_subdirectory(./hydro) #for test case fluid_dynamics ...
add_subdirectory(./jet)  #for test case jet classes …
add_subdirectory(./test)  #for test case ...
add_subdirectory(./reader)  #for test case ...
add_subdirectory(./3rdparty/hydro_from_external_file)  #for test case ...

link_directories(./lib)

add_executable(testJetScape ./src/testJetScape.cc)
add_executable(readerTest ./reader/readerTest.cc)

add_executable(hydroFileTest ./src/hydroFileTest.cc)

target_link_libraries(testJetScape JetScape JetScapeHydro JetScapeTest hydroFromFile) #for test case ...
target_link_libraries(readerTest JetScape JetScapeHydro  JetScapeReader) #for reader test case ...
target_link_libraries(hydroFileTest JetScape JetScapeHydro JetScapeTest hydroFromFile) #for test case ...

if (${HDF5_FOUND})
    set_target_properties (hydroFileTest PROPERTIES COMPILE_FLAGS "-DUSE_HDF5")
    target_link_libraries (hydroFileTest ${_hdf5_libs})
endif()


#Check how to set env variables ...
if (${ROOT_FOUND})
   target_link_libraries(testJetScape ${ROOT_LIBRARIES})
   target_link_libraries(hydroFileTest ${ROOT_LIBRARIES})
endif()

if (${PYTHIA8_FOUND})
    add_executable(writerTest ./src/writerTest.cc)
    add_executable(brickTest ./src/brickTest.cc)
    add_executable(PythiaBrickTest ./src/PythiaBrickTest.cc)
    target_link_libraries(writerTest JetScape JetScapeHydro JetScapeTest) #for writer test case...
    target_link_libraries(brickTest JetScape JetScapeHydro  JetScapeTest) #for test case ...
    target_link_libraries(PythiaBrickTest JetScape JetScapeHydro  JetScapeTest) #for test case ...
    target_link_libraries(testJetScape ${PYTHIA8_LIBRARIES})
    target_link_libraries(brickTest ${PYTHIA8_LIBRARIES})
    target_link_libraries(PythiaBrickTest ${PYTHIA8_LIBRARIES})
    if (${ROOT_FOUND})
      target_link_libraries(brickTest ${ROOT_LIBRARIES})
      target_link_libraries(PythiaBrickTest ${ROOT_LIBRARIES})
      target_link_libraries(writerTest ${ROOT_LIBRARIES})
    endif(${ROOT_FOUND})
    target_link_libraries(hydroFileTest ${PYTHIA8_LIBRARIES})
    target_link_libraries(writerTest ${PYTHIA8_LIBRARIES})
    message( ${PYTHIA8_LIBRARIES} )
endif (${PYTHIA8_FOUND})

# Don't compile tests modules if test is not switched on
if (test)
  # add gtest for the unit test# 
  add_subdirectory("${CMAKE_SOURCE_DIR}/3rdparty/googletest/")
  enable_testing()
  # Include the gtest library. gtest_SOURCE_DIR is available due to
  # 'project(gtest)' above.
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
endif()

if (music)
    message("Building MUSIC ...")
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
    include_directories(./3rdparty/music ./3rdparty/music/src)
    add_subdirectory(./3rdparty/music)
    if (${MPI_FOUND} AND ${GSL_FOUND})
        add_executable(MUSICTest ./src/MUSICTest.cc)
        target_link_libraries(MUSICTest ${MPI_LIBRARIES} ${GSL_LIBRARIES})
        target_link_libraries(testJetScape JetScape JetScapeHydro JetScapeTest hydroFromFile music_lib)
        target_link_libraries(MUSICTest JetScape JetScapeHydro JetScapeTest music_lib)
    endif()
endif (music)

# -----------------------------------------------------------
#### copy jetscape_init.xml to build/, such that one can test xml reading in build/ directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/jetscape_init.xml DESTINATION ${CMAKE_BINARY_DIR}/)

#### same for pythi gun
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/jetscape_init_pythiagun.xml DESTINATION ${CMAKE_BINARY_DIR}/)

#set(LIBS ${LIBS} ${ROOT_LIBRARIES})
#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp")

