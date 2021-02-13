###
cmake_minimum_required (VERSION 2.6) 
project(Foo LANGUAGES CXX)

                             #########
#########################   #Setting#   ##################################
                           #########   
                    
set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_COMPILER "nvc++")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2")

set(SOURCES pricer.cpp 
            hedger.cpp
            gridpricer.cpp
            BSM.сpp
            IRProviderConst.сpp
            )
            
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
            Vanillas.h
            GridOp_S3.h
            GridOp_S3.hpp
            BSM.hpp)
            
set(LIBS IRProviderConst BSM)
            
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

#####################################################################


         #########
######  #Targets#  #######
       #########
set(TARGETS pricer 
            hedger
            gridpricer)

#...
##########################



                             #########
#######################     # loops #     #####################
                           #########

foreach ( LIB ${LIBS} ) 
   add_library(${LIB} ${LIB}.cpp )
endforeach ( LIB )


foreach ( TARGET ${TARGETS} ) 
  add_executable ( ${TARGET} ${TARGET}.cpp )
  target_link_libraries ( ${TARGET} PUBLIC ${LIBS} ) 
endforeach ( TARGET )

#############################################################


                           ###########
########################  # install #  ###################################
#                        ###########                                     #
install(FILES ${SOURCES} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/sources)#
install(FILES ${HEADERS} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include)#
#                                                                        #
##########################################################################


#######Settings flags######
#if(APPLE)

#else()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc=gpu -O2")
#endif()
############


#include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include
#                    ${CMAKE_CURRENT_SOURCE_DIR}/lib
#                    ${CMAKE_CURRENT_SOURCE_DIR}/source
# 

