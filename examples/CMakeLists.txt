#=======================================================================================================================
# Define Demo1 target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-demo1 Demo1.cpp)
target_link_libraries(${PROJECT_NAME}-demo1 PRIVATE common PropertyLib EOSLib FluidLib)

#=======================================================================================================================
# Define Demo2 target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-demo2 Demo2.cpp)
target_link_libraries(${PROJECT_NAME}-demo2 PRIVATE common PropertyLib EOSLib FluidLib)

#=======================================================================================================================
# Define Demo3 target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-demo3 Demo3.cpp)
target_link_libraries(${PROJECT_NAME}-demo3 PRIVATE common PropertyLib EOSLib FluidLib)

#=======================================================================================================================
# Define Demo4 target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-demo4 Demo4.cpp)
target_link_libraries(${PROJECT_NAME}-demo4 PRIVATE common PropertyLib EOSLib FluidLib UnitOps)
