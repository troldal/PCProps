#=======================================================================================================================
# COMMON UTILITIES
#   Define a common utilities library used by HOLD.
#=======================================================================================================================
add_library(pclib INTERFACE IMPORTED GLOBAL)
target_include_directories(pclib SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)


##=======================================================================================================================
## PROJECT FILES
##   List of project source files
##=======================================================================================================================
#set(PROJECT_TEMPLATE_SOURCES
#        PCComponent.cpp
#        )
#
##=======================================================================================================================
## SHARED LIBRARY
##   Define the shared library
##=======================================================================================================================
#add_library(pclib SHARED "")
#target_sources(pclib PRIVATE ${PROJECT_TEMPLATE_SOURCES})
#target_include_directories(pclib
#        PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)     # For export header
#target_link_libraries(pclib PUBLIC dependencies common)
#
#include(GenerateExportHeader)
#generate_export_header(pclib)
##BASE_NAME openxlsx
##EXPORT_FILE_NAME OpenXLSX-Exports.hpp
##EXPORT_MACRO_NAME OPENXLSX_EXPORT
##NO_EXPORT_MACRO_NAME OPENXLSX_HIDDEN)
#
##=======================================================================================================================
## SET TARGET COMPILE OPTIONS
##   ---
##=======================================================================================================================
#
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    target_compile_options(pclib PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNU}>)
#
#elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
#("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") OR
#("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
#    target_compile_options(pclib PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNUCLANG}>)
#
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
#    target_compile_options(pclib PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_MSVC}>)
#    target_compile_options(pclib PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS_MSVC}>)
#endif ()
#
##=======================================================================================================================
## Enable Link-Time Optimization
##=======================================================================================================================
#if (CMAKE_BUILD_TYPE STREQUAL "Release")
#    include(CheckIPOSupported)
#    check_ipo_supported(RESULT result OUTPUT output)
#    if (result)
#        set_property(TARGET pclib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
#    endif ()
#endif ()