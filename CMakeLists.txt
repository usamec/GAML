cmake_minimum_required (VERSION 2.6)
project (Ass)

FIND_PACKAGE( Boost 1.46 COMPONENTS serialization REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

#SET(GCC_COVERAGE_COMPILE_FLAGS "--coverage")
#SET(GCC_COVERAGE_LINK_FLAGS    "--coverage")

list( APPEND CMAKE_CXX_FLAGS "-std=c++0x -g -O2 ${CMAKE_CXX_FLAGS}")

add_library(graph graph.cc)

add_library(input_output input_output.cc)
target_link_libraries(input_output graph)

add_library(moves moves.cc)
target_link_libraries(moves graph)

add_library(graph_from_assembly graph_from_assembly.cc)
target_link_libraries(graph_from_assembly graph)

add_executable(gaml gaml.cc)
target_link_libraries(gaml graph input_output moves graph_from_assembly ${Boost_LIBRARIES})

add_executable(testret testrep.cc)
