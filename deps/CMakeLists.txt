#=======================================================================================================================
# EXTERNAL LIBRARIES
#   Define external libraries used by PCProps. The libraries are header-only, so
#   INTERFACE libraries should be defined.
#=======================================================================================================================
add_library(dependencies INTERFACE)
target_include_directories(dependencies SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
