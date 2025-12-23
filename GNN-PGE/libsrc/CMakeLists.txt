# add graph lib source code
file(GLOB GRAPH_FILE ${PROJECT_SOURCE_DIR}/libsrc/graph/*.cpp)
message(STATUS "GRAPH_SOURCES: " ${GRAPH_FILE})
add_library(graph OBJECT  ${GRAPH_FILE})

# add blockfile lib source code
file(GLOB BLOCKFILE_FILE ${PROJECT_SOURCE_DIR}/libsrc/blockfile/*.cpp)
message(STATUS "BLOCKFILE_SOURCES: " ${BLOCKFILE_FILE})
add_library(blockfile OBJECT  ${BLOCKFILE_FILE})

# add func lib source code
file(GLOB FUNC_FILE ${PROJECT_SOURCE_DIR}/libsrc/func/*.cpp)
message(STATUS "FUNC_SOURCES: " ${FUNC_FILE})
add_library(func OBJECT  ${FUNC_FILE})

# add heap lib source code
file(GLOB HEAP_FILE ${PROJECT_SOURCE_DIR}/libsrc/heap/*.cpp)
message(STATUS "HEAP_SOURCES: " ${HEAP_FILE})
add_library(heap OBJECT  ${HEAP_FILE})

# add linlist source code
file(GLOB LINLIST_FILE ${PROJECT_SOURCE_DIR}/libsrc/linlist/*.cpp)
message(STATUS "LINLIST_SOURCES: " ${LINLIST_FILE})
add_library(linlist OBJECT  ${LINLIST_FILE})

# add utility source code
file(GLOB UTILITY_FILE ${PROJECT_SOURCE_DIR}/libsrc/utility/*.cpp
${PROJECT_SOURCE_DIR}/libsrc/utility/han/*.cpp 
${PROJECT_SOURCE_DIR}/libsrc/utility/han/utils/*.cpp)
message(STATUS "UTILITY_SOURCES: " ${UTILITY_FILE})
add_library(utility OBJECT  ${UTILITY_FILE})

# add rtree source code
file(GLOB RTREE_FILE ${PROJECT_SOURCE_DIR}/libsrc/rtree/*.cpp)
message(STATUS "RTREE_SOURCES: " ${RTREE_FILE})
add_library(rtree OBJECT ${RTREE_FILE})