#=======================================================================================================================
# PROJECT FILES
#   List of project source files
#=======================================================================================================================
set(PROPERTIESLIB_SOURCES
        PhaseProperties.cpp
        FluidProperties.cpp
        )

#=======================================================================================================================
# SHARED LIBRARY
#   Define the shared library
#=======================================================================================================================
add_library(PropertiesLib STATIC "")
target_sources(PropertiesLib PRIVATE ${PROPERTIESLIB_SOURCES})
target_include_directories(PropertiesLib
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)     # For export header
target_link_libraries(PropertiesLib PUBLIC dependencies)

#include(GenerateExportHeader)
#generate_export_header(FluidLib)
#BASE_NAME openxlsx
#EXPORT_FILE_NAME OpenXLSX-Exports.hpp
#EXPORT_MACRO_NAME OPENXLSX_EXPORT
#NO_EXPORT_MACRO_NAME OPENXLSX_HIDDEN)

#=======================================================================================================================
# SET TARGET COMPILE OPTIONS
#   ---
#=======================================================================================================================

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(PropertiesLib PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNU}>)

elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    target_compile_options(PropertiesLib PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNUCLANG}>)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    target_compile_options(PropertiesLib PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_MSVC}>)
    target_compile_options(PropertiesLib PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS_MSVC}>)
endif ()

#=======================================================================================================================
# Enable Link-Time Optimization
#=======================================================================================================================
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output)
    if (result)
        set_property(TARGET PropertiesLib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif ()
endif ()