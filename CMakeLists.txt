###
cmake_minimum_required (VERSION 2.6) 
project(Foo)

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2")
#else()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
#endif()
############


#######Targets#########
##set(TARGETS pricer hedger)
#...
#######################
add_executable(pricer main.cpp )

add_executable(hedger mainhedge.cpp )

add_library(IRProvider IRProvider.cpp)

target_link_libraries(pricer PUBLIC IRProvider)
target_link_libraries(hedger PUBLIC IRProvider)

install(TARGETS IRProvider DESTINATION lib)
install(TARGETS pricer hedger DESTINATION bin)
install(FILES *.h *.hpp DESTINATION include)
install(FILES *.cpp DESTINATION source)



