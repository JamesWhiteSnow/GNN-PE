cmake_minimum_required(VERSION 3.0) 

project(rtree) 

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++17 -O3 -g -march=native -pthread -fopenmp") 

include_directories(${PROJECT_SOURCE_DIR}/include ) 

add_subdirectory(${PROJECT_SOURCE_DIR}/libsrc)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

