############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cap_test CXX)

add_definitions(-pg)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

include_directories(${CMAKE_SOURCE_DIR}/debruijn)
include_directories(${CMAKE_SOURCE_DIR}/cap)

add_executable(cap_test
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/kmer_coverage_model.cpp
 test.cpp)

target_link_libraries(cap_test input mph_index otkpp ${COMMON_LIBRARIES})

