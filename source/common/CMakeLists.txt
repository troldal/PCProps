#=======================================================================================================================
# COMMON UTILITIES
#   Define a common utilities library used by HOLD.
#=======================================================================================================================
add_library(common INTERFACE IMPORTED GLOBAL)
target_include_directories(common SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)