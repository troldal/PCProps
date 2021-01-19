#=======================================================================================================================
# COMMON UTILITIES
#   Define a common utilities library used by HOLD.
#=======================================================================================================================
add_library(pclib INTERFACE)
target_include_directories(pclib SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
