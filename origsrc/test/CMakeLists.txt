##############
# Unit Tests
##############

###### Begin: Helper macro to add unittests
macro(js_add_test name depends)
   add_test(NAME ${name} COMMAND ${ARGN})
   add_custom_target(run_${name}
      COMMAND ${ARGN} ${UNITTEST_ARGS}
      DEPENDS ${depends}
      COMMENT "Executing test ${name}"
      VERBATIM
      )
endmacro()

macro(js_add_exe name)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} gtest gtest_main)
  target_link_libraries(${name} JetScape)
endmacro()

macro(js_add_unittest name)
   js_add_exe(${name})
   js_add_test(${name} ${name} ${name})
endmacro()
###### End: Helper macro to add unittests


###### Now it is easy to add unittest for test_name.cc by js_add_unittest(test_name)
js_add_unittest(TrentoInitial)
js_add_unittest(TrentoInitialFromXML)
