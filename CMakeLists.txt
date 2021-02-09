###
cmake_minimum_required (VERSION 2.6) 
project(foo)

#######Settings flags######
if(APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()
############

########Setting#######
set(CMAKE_CXX_COMPILER /bin/clang++)
set(HEADERS Diffusion.h 
            IRProvider.h 
            IRProviderConst.h 
            MonteCarlo.h 
            MonteCarlo.hpp)
set(SOURCES main.cpp 
            IRProvider.cpp)
######################

#######Targets#########

#...
#######################

add_executable(foo $(HEADERS) $(SOURCES))


