##############
# Unit Tests
##############
add_executable(runUnitTests linear_interpolation.cc)

# Standard linking to gtest stuff.
target_link_libraries(runUnitTests gtest gtest_main)

# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable runUnitTests to see those specific tests.
add_test(NAME test1 COMMAND runUnitTests)

# You can also omit NAME and COMMAND. The second argument could be some other
# test executable.
add_test(NAME test2 COMMAND runUnitTests)
