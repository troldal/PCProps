#=======================================================================================================================
# Interface Modules
#   Define a common interface modules
#=======================================================================================================================
add_library(interfaces INTERFACE IMPORTED GLOBAL)
target_include_directories(interfaces SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)