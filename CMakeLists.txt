cmake_minimum_required (VERSION 2.6)
project (Ass)

FIND_PACKAGE( Boost 1.46 COMPONENTS serialization REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

list( APPEND CMAKE_CXX_FLAGS "-std=c++0x -g -O2 ${CMAKE_CXX_FLAGS}")

add_library(graph graph.cc)

add_library(input_output input_output.cc)
target_link_libraries(input_output graph)

add_library(moves moves.cc)
target_link_libraries(moves graph)

#add_executable(ass ass.cc)
#target_link_libraries(ass graph ${Boost_LIBRARIES})

#add_executable(ass2 ass2.cc)
#target_link_libraries(ass2 graph ${Boost_LIBRARIES})
#
#add_executable(ass454 ass454.cc)
#target_link_libraries(ass454 graph ${Boost_LIBRARIES})
#
#add_executable(assshort assshort.cc)
#target_link_libraries(assshort graph ${Boost_LIBRARIES})
#
#add_executable(postprocess postprocess.cc)
#target_link_libraries(postprocess graph ${Boost_LIBRARIES})
#
#add_executable(postprocess2 postprocess2.cc)
#target_link_libraries(postprocess2 graph ${Boost_LIBRARIES})
#
#add_executable(postprocessr postprocessr.cc)
#target_link_libraries(postprocessr graph ${Boost_LIBRARIES})
#
#add_executable(postprocess454 postprocess454.cc)
#target_link_libraries(postprocess454 graph ${Boost_LIBRARIES})
#
#add_executable(testp testp.cc)
#target_link_libraries(testp graph ${Boost_LIBRARIES})
#
#add_executable(calcprob calcprob.cc)
#target_link_libraries(calcprob graph ${Boost_LIBRARIES})
#
#add_executable(ssc ssc.cc)
#target_link_libraries(ssc graph ${Boost_LIBRARIES})
#
#add_executable(prefilter prefilter.cc)
#target_link_libraries(prefilter ${Boost_LIBRARIES})

add_executable(lama lama.cc)
target_link_libraries(lama graph input_output moves ${Boost_LIBRARIES})

