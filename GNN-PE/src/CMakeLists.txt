aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)

add_executable(main ${SRC_LIST})

 
target_link_libraries(main 
func 
graph 
blockfile 
heap 
linlist
rtree
utility
) 