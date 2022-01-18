add_library(sciplot INTERFACE IMPORTED GLOBAL)
target_include_directories(sciplot SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/plot)


#=======================================================================================================================
# Define Demo1 target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-demo1 Demo1.cpp)
target_link_libraries(${PROJECT_NAME}-demo1 PRIVATE common PureComponentLib PropertiesLib PropertyLib PropertyPackage FluidLib sciplot)

#=======================================================================================================================
# Define Demo2 target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-demo2 Demo2.cpp)
target_link_libraries(${PROJECT_NAME}-demo2 PRIVATE common PureComponentLib PropertiesLib PropertyLib PropertyPackage FluidLib sciplot UnitOps)

