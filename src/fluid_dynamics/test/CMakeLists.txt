##############
# Unit Tests
##############
add_executable(test_linear_interpolation1d linear_interpolation.cc)

# Standard linking to gtest stuff.
target_link_libraries(test_linear_interpolation1d gtest gtest_main)

# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable runUnitTests to see those specific tests.
add_test(NAME test1 COMMAND test_linear_interpolation1d)
