# Compile everything except the main source file into a static lib to be linked
# to both the main executable and the tests.
add_library(${LIBRARY_NAME} STATIC
  collider.cxx
  event.cxx
  hdf5_utils.cxx
  nucleon.cxx
  nucleus.cxx
  output.cxx
  random.cxx
)
set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")

# Compile the actual executable.
set(MAIN trento.cxx)
set_source_files_properties(${MAIN} PROPERTIES
  COMPILE_DEFINITIONS "TRENTO_VERSION_STRING=\"${PROJECT_VERSION}\"")
add_executable(${PROJECT_NAME} ${MAIN})
target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME} ${Boost_LIBRARIES} ${HDF5_LIBRARIES} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
