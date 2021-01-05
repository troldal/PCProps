#=======================================================================================================================
# Define Demo1 target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-demo1 Demo1.cpp)
target_link_libraries(${PROJECT_NAME}-demo1 PRIVATE ${PROJECT_NAME}-library-static)
target_compile_definitions(${PROJECT_NAME}-demo1 PRIVATE OPENXLSX_STATIC_DEFINE)

#=======================================================================================================================
# Define Demo2 target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-demo2 Demo2.cpp)
target_link_libraries(${PROJECT_NAME}-demo2 PRIVATE ${PROJECT_NAME}-library-static)
target_compile_definitions(${PROJECT_NAME}-demo2 PRIVATE OPENXLSX_STATIC_DEFINE)

