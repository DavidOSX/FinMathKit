###
cmake_minimum_required (VERSION 2.6) 
project(Foo)

########Setting#######
set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_COMPILER "nvc++")
set(SOURCES pricer.cpp 
            IRProvider.cpp 
            hedger.cpp)
            
set(HEADERS Diffusion.h 
            IRProvider.h
            IRProviderConst.h
            MCOptionHedger.h
            MCOptionHedger.hpp
            MCOptionPricer.h
            MCOptionPricer.hpp
            MonteCarlo.h
            MonteCarlo.hpp
            Options.h
            Time.h
            Vanillas.h)
            
set(LIBS IRProvider)
            
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

######################


#######Settings flags######
#if(APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2")
#else()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc=gpu -O2")
#endif()
############


#######Targets#########
set(TARGETS pricer 
            hedger)

#...
#######################
#include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include
#                    ${CMAKE_CURRENT_SOURCE_DIR}/lib
#                    ${CMAKE_CURRENT_SOURCE_DIR}/source
#                    )


add_library(${LIBS} IRProvider.cpp)

foreach ( TARGET ${TARGETS} ) 
  add_executable ( ${TARGET} ${TARGET}.cpp )
  target_link_libraries ( ${TARGET} PUBLIC ${LIBS} ) 
endforeach ( TARGET )

###########################
install(FILES ${SOURCES} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/sources)
install(FILES ${HEADERS} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include)
###########################



