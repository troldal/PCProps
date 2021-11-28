#=======================================================================================================================
# PROJECT FILES
#   List of project source files
#=======================================================================================================================
set(EOSLIB_SOURCES
        numerics.cpp
        )

#=======================================================================================================================
# SHARED LIBRARY
#   Define the shared library
#=======================================================================================================================
set(PCPROPS_USE_GSL ON)

add_library(numerics STATIC "")
target_sources(numerics PRIVATE ${EOSLIB_SOURCES})
target_include_directories(numerics
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)     # For export header

if (${PCPROPS_USE_GSL})
    find_package(GSL REQUIRED)
    target_link_libraries(numerics PUBLIC GSL::gsl)
    add_compile_definitions(PCPROPS_USE_GSL)
endif()

#include(GenerateExportHeader)
#generate_export_header(EOSLib)
#BASE_NAME openxlsx
#EXPORT_FILE_NAME OpenXLSX-Exports.hpp
#EXPORT_MACRO_NAME OPENXLSX_EXPORT
#NO_EXPORT_MACRO_NAME OPENXLSX_HIDDEN)

#=======================================================================================================================
# SET TARGET COMPILE OPTIONS
#   ---
#=======================================================================================================================

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(numerics PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNU}>)

elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    target_compile_options(numerics PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNUCLANG}>)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    target_compile_options(numerics PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_MSVC}>)
    target_compile_options(numerics PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS_MSVC}>)
endif ()

#=======================================================================================================================
# Enable Link-Time Optimization
#=======================================================================================================================
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output)
    if (result)
        set_property(TARGET numerics PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif ()
endif ()