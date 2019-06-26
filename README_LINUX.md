In order to install and compile the JETSCAPE framework on a Linux machine, one need to follow the following steps:

1. Clone the repository from Github:
     - git clone https://github.com/amajumder/JETSCAPE-COMP.git
     - Enter your Github credentials 

2. Install and configure HepMC (Version 3.0.0 or higher)
     - Create a folder for installing HepMC
     - Go to the created folder
     - wget http://hepmc.web.cern.ch/hepmc/releases/hepmc3.0.0.tgz
     - tar -xvzf hepmc3.0.0.tgz
     - Go to the extracted folder 
          - cd hepmc3.0.0/cmake
     - Run cmake command
          - cmake ..
     - Run make command (the user has to be sudoer)
          - make all install
     - The default path to install HepMC is /usr/local/lib/ 

3. Change Cmake configuration for HepMC (If Cmake could not find HepMC)
     - Set "HEPMC_DIR" as the root directory of HepMC
     - Go to ../framework/Modules
     - Open FindHEPMC.cmake
          - vim FindHEPMC.cmake
     - Make sure that cmake looks for the correct "include" and "lib" directories of HepMC
          - By default is "${HEPMC_DIR}/include/"
          - By default is "${HEPMC_DIR}/lib"
          - Replace "${HEPMC_DIR}" with "ENV{HEPMC_DIR}" if you are not root
     - The library file that cmake must find is "libHepMC.so"


4. Install Pythia8
     - Create a folder for installing Pythia8
     - Got to the created folder
     - wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8226.tgz
     - tar -xvzf pythia8226.tgz
     - Go to the extracted folder 
          - cd pythia8226
     - Run make command
          - make

5. Configure the address of Pythia8 in the ../framework/setup.sh (If cmake could not find Pythia)
     - vim setup.sh
     - Change "PYTHIAINSTALLDIR" to the installing folder of Pythia8
     - Update the "PYTHIA8DIR" and "PYTHIA8_ROOT_DIR" address
          - should be "${PYTHIAINSTALLDIR}/pythia8226" by default
     - ./setup.sh
     - Make sure the variable are set into the session
     - Use source if they are not set

6. Install Boost libraries (Version 1.5 or higher)
     - wget https://dl.bintray.com/boostorg/release/1.64.0/source/boost_1_64_0.tar.gz 
     - tar -xvzf boost_1_64_0.tar.gz
     - Go to the directory tools/build/.
     - Run bootstrap.sh 
          - ./bootstrap.sh 
     - Run b2 --prefix=PREFIX where PREFIX is the directory where you want Boost. Build to be installed
          - ./b2 install --prefix=PREFIX
     - Set "BOOST_ROOT" variable to the root directory of boost
     - export BOOST_ROOT=<root_directory> 

7. Install HDF5 C++ library
     - wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.18/src/hdf5-1.8.18.tar
     - tar -xvvf hdf5-1.8.18.tar
     - ./configure
     - make

8. Create a build folder in ../src/framework
     - mkdir build
     - Got to the build folder
          - cd build
     - Make sure that in CMakeLists.text there is a condition on UNIX not LINUX
          - if(UNIX) 
               message( STATUS "Linux : " ${CMAKE_HOST_SYSTEM})
               set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
            endif()
     - Make sure $LD_LIBRARY_PATH is set to the address of dynamic library files
          - echo $LD_LIBRARY_PATH
          - If there is nothing to dispay
               - export LD_LIBRARY_PATH=<dynamic_lib_dir>
          - If the address is not included in the variable
               - LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<dynamic_lib_dir>
     - Run cmake
          - cmake ..
          - If cmake cannot find HDF5 library, set "-DCMAKE_LIBRARY_PATH" and "-DCMAKE_INCLUDE_PATH" flags when you run cmake
               - cmake -DCMAKE_LIBRARY_PATH=<HDF5_Include_Path> -DCMAKE_INCLUDE_PATH=<HDF5_Lib_Path> ..
     - run make
          - make

