##############
# Unit Tests
##############
add_executable(test_jet jet_class.cc)

# Standard linking to gtest stuff.
target_link_libraries(test_jet gtest gtest_main)
target_link_libraries(test_jet jetclass_lib)

# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable test_jet to see those specific tests.
add_test(NAME test1 COMMAND test_jet)
