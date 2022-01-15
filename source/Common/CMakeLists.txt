#=======================================================================================================================
# COMMON UTILITIES
#   Define a Common utilities library used by HOLD.
#=======================================================================================================================
add_library(common INTERFACE IMPORTED GLOBAL common/Globals.hpp common/PCPropsException.hpp)
target_include_directories(common SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

