##############
# Unit Tests
##############
add_executable(test_linear_interpolation1d linear_interpolation.cc)
add_executable(test_fluid_dynamics fluid_dynamics.cc)

# Standard linking to gtest stuff.
target_link_libraries(test_linear_interpolation1d gtest gtest_main)
target_link_libraries(test_fluid_dynamics gtest gtest_main)
target_link_libraries(test_fluid_dynamics hydro_lib)

# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable runUnitTests to see those specific tests.
add_test(NAME test1 COMMAND test_linear_interpolation1d)
add_test(NAME test2 COMMAND test_fluid_dynamics)
