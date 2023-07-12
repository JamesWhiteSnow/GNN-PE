cmake_minimum_required(VERSION 3.0) 

project(rtree) 

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -o3 -mavx2 -w --std=c++11") 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g -march=native -pthread") 

# add_compile_definitions(UNIX)

include_directories(${PROJECT_SOURCE_DIR}/include ) 

add_subdirectory(${PROJECT_SOURCE_DIR}/libsrc)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

