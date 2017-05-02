##############
# Unit Tests
##############
add_executable(test_jetscape_ini JetScapeInitial.cc)

# Standard linking to gtest stuff.
target_link_libraries(test_jetscape_ini gtest gtest_main)
# trento asks for Boost and HDF5 libraries
target_link_libraries(test_jetscape_ini JetScape)


# trento asks for Boost and HDF5 libraries
#target_link_libraries(test_jetscape_ini libtrento ${Boost_LIBRARIES} ${HDF5_LIBRARIES})
# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable runUnitTests to see those specific tests.
add_test(NAME test1 COMMAND test_jetscape_ini)
