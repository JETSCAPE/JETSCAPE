find_program(DOXYGEN doxygen)
find_program(SPHINX sphinx-build)

if(DOXYGEN AND SPHINX)
  configure_file(Doxyfile.in Doxyfile)
  configure_file(conf.py.in conf.py)
  add_custom_target(doc
    COMMAND ${DOXYGEN}
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMAND ${SPHINX}
      -b html
      -c ${CMAKE_CURRENT_BINARY_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}/html
  )
endif()
