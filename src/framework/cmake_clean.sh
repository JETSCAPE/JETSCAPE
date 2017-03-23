#bin/sh

rm CMakeCache.txt
rm -rf CMakeFiles
rm cmake_install.cmake
rm Makefile

cd ./src

rm -rf CMakeFiles
rm cmake_install.cmake
rm Makefile

cd ../test

rm -rf CMakeFiles
rm cmake_install.cmake
rm Makefile

cd ../hydro

rm -rf CMakeFiles
rm cmake_install.cmake
rm Makefile

cd ../jet

rm -rf CMakeFiles
rm cmake_install.cmake
rm Makefile
