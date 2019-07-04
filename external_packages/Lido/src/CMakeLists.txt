# compile everything except the main source file into a statis lib to be linked 
# to the main executable
add_library(${LIBRARY_NAME} STATIC
TableBase.cpp
StochasticBase.cpp
Xsection.cpp
Rate.cpp
matrix_elements_1-2.cpp
matrix_elements_2-2.cpp
matrix_elements_2-3.cpp
hcubature.cpp
random.cpp
simpleLogger.cpp
stat.cpp
approx_functions.cpp
workflow.cpp
Langevin.cpp
predefine.cpp
lorentz.cpp
Medium_Reader.cpp
	)

set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")


# compile the actual executable
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME} ${PYTHIA8_LIBRARIES} ${GSL_LIBRARIES} ${GSLCALAS_LIBRARIES} ${HDF5_LIBRARIES} ${Boost_LIBRARIES} -pthread -lpthread)

# install executable
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
# also install the generated c++ shared library
install(TARGETS ${LIBRARY_NAME} DESTINATION lib/${LIBRARY_NAME})
# install the library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/${LIBRARY_NAME})
