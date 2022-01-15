#=======================================================================================================================
# COMMON UTILITIES
#   Define a Common utilities library used by HOLD.
#=======================================================================================================================
add_library(PropertyLib INTERFACE)
target_include_directories(PropertyLib SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
