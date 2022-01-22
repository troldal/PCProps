#=======================================================================================================================
# PROJECT FILES
#   List of project source files
#=======================================================================================================================
set(UNITOPS_SOURCES
#        Pipe/Pipe.cpp
#        Pipe/PipeSegmentLiquid.cpp

        Pump/CentrifugalPump.cpp
        Stream/Stream.cpp
        )

#=======================================================================================================================
# SHARED LIBRARY
#   Define the shared library
#=======================================================================================================================
add_library(UnitOps SHARED "")
target_sources(UnitOps PRIVATE ${UNITOPS_SOURCES})
target_include_directories(UnitOps
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)     # For export header
target_link_libraries(UnitOps PUBLIC dependencies common interfaces PropertiesLib)

include(GenerateExportHeader)
generate_export_header(UnitOps)
#BASE_NAME openxlsx
#EXPORT_FILE_NAME OpenXLSX-Exports.hpp
#EXPORT_MACRO_NAME OPENXLSX_EXPORT
#NO_EXPORT_MACRO_NAME OPENXLSX_HIDDEN)

#=======================================================================================================================
# SET TARGET COMPILE OPTIONS
#   ---
#=======================================================================================================================

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(UnitOps PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNU}>)

elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    target_compile_options(UnitOps PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNUCLANG}>)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    target_compile_options(UnitOps PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_MSVC}>)
    target_compile_options(UnitOps PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS_MSVC}>)
endif ()

#=======================================================================================================================
# Enable Link-Time Optimization
#=======================================================================================================================
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output)
    if (result)
        set_property(TARGET UnitOps PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif ()
endif ()