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

message ("Include Directories ...")

include_directories(./src ./test ./hydro ./jet ./gtl/include) #for test case ...

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
     message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})	
     message(STATUS "LD_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})
endif()   

# w/o not compiling library automatically
add_subdirectory(./gtl)
add_subdirectory(./src)
add_subdirectory(./hydro) #for test case fluid_dynamics ...
add_subdirectory(./jet)  #for test case jet classes …
add_subdirectory(./test)  #for test case ...

link_directories(./lib)

add_executable(testJetScape ./src/testJetScape.cc)
add_executable(brickTest ./src/brickTest.cc)

target_link_libraries(testJetScape JetScape JetScapeHydro JetScapeTest) #for test case ...
target_link_libraries(brickTest JetScape JetScapeHydro  JetScapeTest) #for test case ...

#Check how to set env variables ...
if (${ROOT_FOUND})
   target_link_libraries(testJetScape ${ROOT_LIBRARIES})
   target_link_libraries(brickTest ${ROOT_LIBRARIES})
endif()

if (${PYTHIA8_FOUND})
   target_link_libraries(testJetScape ${PYTHIA8_LIBRARIES})
   target_link_libraries(brickTest ${PYTHIA8_LIBRARIES})
endif()

#add_executable(testJS ./src/testJetScape.cc)
#target_link_libraries(testJS JetScape)

# -----------------------------------------------------------
    
#set(LIBS ${LIBS} ${ROOT_LIBRARIES})
#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp")
