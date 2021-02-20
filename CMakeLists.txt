###
cmake_minimum_required (VERSION 3.0) 
project(Foo LANGUAGES CXX)

                             #########
#########################   #Setting#   ##################################
                           #########   
                    
set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_COMPILER "nvc++")
#set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-10")

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

                   ###########
################# # testing # #################
                 ###########
                 
enable_testing()

function(do_test1 TARGET ARG1 ARG2 ARG3 ARG4 ARG6 ARG7 ARG8 ARG9 )
  add_test(pricer:Test1 ${TARGET} ${ARG1} ${ARG2} ${ARG3} ${ARG4} "Call" ${ARG6} ${ARG7} ${ARG8} ${ARG9} "GBM" 0 )
  
  add_test(pricer:Test2 ${TARGET} ${ARG1} ${ARG2} ${ARG3} ${ARG4} "Put" ${ARG6} ${ARG7} ${ARG8} ${ARG9} "GBM" 0 )
  
  add_test(pricer:Test3 ${TARGET} ${ARG1} ${ARG2} ${ARG3} ${ARG4} "Call" ${ARG6} ${ARG7} ${ARG8} ${ARG9} "CEV" 0 )
  
  add_test(pricer:Test4 ${TARGET} ${ARG1} ${ARG2} ${ARG3} ${ARG4} "Put" ${ARG6} ${ARG7} ${ARG8} ${ARG9} "CEV" 0 )
  #set_tests_properties(Test1 ...)
endfunction(do_test1)

function(do_test2 TARGET ARG1 ARG2 ARG3 ARG4 ARG6 ARG7 ARG8 ARG9 ARG10 )
  add_test(hedger:Test1 ${TARGET} ${ARG1} ${ARG2} ${ARG3} ${ARG4} "Call" ${ARG6} ${ARG7} ${ARG8} ${ARG9} ${ARG10} "GBM" 0 )
  
  add_test(hedger:Test2 ${TARGET} ${ARG1} ${ARG2} ${ARG3} ${ARG4} "Put" ${ARG6} ${ARG7} ${ARG8} ${ARG9} ${ARG10} "GBM" 0 )
  
  #add_test(hedger:Test3 ${TARGET} ${ARG1} ${ARG2} ${ARG3} ${ARG4} "Call" ${ARG6} ${ARG7} ${ARG8} ${ARG9} ${ARG10} "CEV" 0 )
  #set_tests_properties(Test2 ...) 
endfunction(do_test2)

function(do_test3 TARGET ARG1 ARG2 ARG3 ARG5 ARG6 ARG7 ARG8 )
  add_test(gridpricer:Test1 ${TARGET} ${ARG1} ${ARG2} ${ARG3} "Call" ${ARG5} ${ARG6} ${ARG7} ${ARG8} "GBM" 0 )
  
  add_test(gridpricer:Test2 ${TARGET} ${ARG1} ${ARG2} ${ARG3} "Put" ${ARG5} ${ARG6} ${ARG7} ${ARG8} "GBM" 0 )
  
  add_test(gridpricer:Test3 ${TARGET} ${ARG1} ${ARG2} ${ARG3} "Call" ${ARG5} ${ARG6} ${ARG7} ${ARG8} "GBM" 1 )
  
  add_test(gridpricer:Test4 ${TARGET} ${ARG1} ${ARG2} ${ARG3} "Put" ${ARG5} ${ARG6} ${ARG7} ${ARG8} "GBM" 1 )
  
  add_test(gridpricer:Test5 ${TARGET} ${ARG1} ${ARG2} ${ARG3} "Put" ${ARG5} ${ARG6} ${ARG7} ${ARG8} "CEV" 0 )
  
  add_test(gridpricer:Test6 ${TARGET} ${ARG1} ${ARG2} ${ARG3} "Put" ${ARG5} ${ARG6} ${ARG7} ${ARG8} "CEV" 1 )
  #set_tests_properties(Test3 ...)
endfunction(do_test3)

do_test1(${CMAKE_CURRENT_SOURCE_DIR}/bin/pricer "DATA.txt" 0 0.25 100 100 90 30 2048)

do_test2(${CMAKE_CURRENT_SOURCE_DIR}/bin/hedger "DATA.txt" 0 0.25 100 100 90 0.01 30 2048)

do_test3(${CMAKE_CURRENT_SOURCE_DIR}/bin/gridpricer "DATA.txt" 0.25 100 100 90 1024 5)
                 
##############################################

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

