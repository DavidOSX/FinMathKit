###
cmake_minimum_required (VERSION 2.6) 
project(foo)

#######Settings flags######
if(APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()
############

#######Targets#########
set(TARGETS foo)
#...
#######################

add_executable(foo main.cpp IRProvider.cpp)

#install( TARGETS ${TARGETS}
         RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin/ )
