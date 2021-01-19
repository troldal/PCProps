#=======================================================================================================================
# COMMON UTILITIES
#   Define a common utilities library used by HOLD.
#=======================================================================================================================
add_library(PropertyLib INTERFACE)
target_include_directories(PropertyLib SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
