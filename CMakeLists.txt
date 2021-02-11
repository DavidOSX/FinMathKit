###
cmake_minimum_required (VERSION 2.6) 
project(foo)

########Setting#######
set(CMAKE_CXX_COMPILER "clang++")
#set(HEADERS Diffusion.h 
#            IRProvider.h 
#            IRProviderConst.h 
#            MonteCarlo.h 
#            MonteCarlo.hpp)
#set(SOURCES main.cpp 
#            IRProvider.cpp)
######################


#######Settings flags######
#if(APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
#else()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
#endif()
############


#######Targets#########

#...
#######################

add_executable(foo main.cpp IRProvider.cpp)


