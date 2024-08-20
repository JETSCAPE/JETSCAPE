# ##############################################################################
# Unit Tests
# ##############################################################################

# tests: Helper macro (from smash) that adds a run_<name> target
macro(js_add_test name depends)
  add_test(NAME ${name} COMMAND ${ARGN})
  add_custom_target(
    run_${name}
    COMMAND ${ARGN} ${UNITTEST_ARGS}
    DEPENDS ${depends}
    COMMENT "Executing test ${name}"
    VERBATIM)
endmacro()

macro(js_add_exe name)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} JetScape)
  target_link_libraries(${name} gtest gtest_main)
  target_link_libraries(${name} gtest gtest_main)
endmacro()

# this is the main func called by developer to add unittests
macro(add_unittest name)
  js_add_exe(${name})
  js_add_test(${name} ${name} ${name})
endmacro()

# add your unittests after this line

add_unittest(linear_interpolation)
add_unittest(fluid_dynamics)
add_unittest(causal_liquifier)
add_unittest(LiquifierBase)
