#=======================================================================================================================
# COMPILER FLAGS
#   Set compiler debug flags for GCC, Clang and MSVC.
#=======================================================================================================================
list(APPEND DEBUG_FLAGS_GNU
        "-Wmisleading-indentation"
        "-Wduplicated-cond"
        "-Wduplicated-branches"
        "-Wlogical-op"
        "-Wnull-dereference"
        "-Wuseless-cast")
list(APPEND DEBUG_FLAGS_GNUCLANG
        "-Wall"
        "-Wextra"
        "-Wshadow"
        "-Wnon-virtual-dtor"
        "-Wold-style-cast"
        "-Wcast-align"
        "-Wunused"
        "-Woverloaded-virtual"
        "-Wpedantic"
        "-Wconversion"
        "-Wsign-conversion"
        "-Wdouble-promotion"
        "-Wformat=2"
        "-Weffc++"
        "-Wno-unknown-pragmas")
list(APPEND DEBUG_FLAGS_MSVC
        "/permissive"
        "/W4"
        "/w14242"
        "/w14254"
        "/w14263"
        "/w14265"
        "/w14287"
        "/we4289"
        "/w14296"
        "/w14311"
        "/w14545"
        "/w14546"
        "/w14547"
        "/w14549"
        "/w14555"
        "/w14619"
        "/w14640"
        "/w14826"
        "/w14905"
        "/w14906"
        "/w14928"
        "/wd4251"
        "/wd4275")

list(APPEND RELEASE_FLAGS_MSVC
        "/wd4251"
        "/wd4275")


add_subdirectory(common)
add_subdirectory(eos)
add_subdirectory(pcomp)

##=======================================================================================================================
## PROJECT FILES
##   List of project source files
##=======================================================================================================================
#set(PROJECT_TEMPLATE_SOURCES
#        pcomp/PCComponent.cpp
#        pcomp/ConstantData/CDJoback.cpp
#        eos/PengRobinson/PengRobinson.cpp
#        pcomp/UnitOps/Pipe.cpp
#        )
#
#
##=======================================================================================================================
## STATIC LIBRARY
##   Define the static library
##=======================================================================================================================
#add_library(${PROJECT_NAME}-library-static STATIC "")
#target_sources(${PROJECT_NAME}-library-static PRIVATE ${PROJECT_TEMPLATE_SOURCES})
#target_include_directories(${PROJECT_NAME}-library-static
#        PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/library>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)     # For export header
#target_link_libraries(${PROJECT_NAME}-library-static PUBLIC dependencies common)
#target_compile_definitions(${PROJECT_NAME}-library-static PUBLIC OPENXLSX_STATIC_DEFINE)
#
##=======================================================================================================================
## SHARED LIBRARY
##   Define the shared library
##=======================================================================================================================
#add_library(${PROJECT_NAME}-library-shared SHARED "")
#target_sources(${PROJECT_NAME}-library-shared PRIVATE ${PROJECT_TEMPLATE_SOURCES})
#target_include_directories(${PROJECT_NAME}-library-shared
#        PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/library>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)     # For export header
#target_link_libraries(${PROJECT_NAME}-library-shared PUBLIC dependencies common)
#
##=======================================================================================================================
## GENERATE EXPORT HEADER
##   Generate preprocessor definitions for setting symbol visibility.
##=======================================================================================================================
#include(GenerateExportHeader)
#generate_export_header(${PROJECT_NAME}-library-shared)
#        #BASE_NAME openxlsx
#        #EXPORT_FILE_NAME OpenXLSX-Exports.hpp
#        #EXPORT_MACRO_NAME OPENXLSX_EXPORT
#        #NO_EXPORT_MACRO_NAME OPENXLSX_HIDDEN)
#
##=======================================================================================================================
## SET TARGET COMPILE OPTIONS
##   ---
##=======================================================================================================================
#
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    target_compile_options(${PROJECT_NAME}-library-shared PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNU}>)
#    target_compile_options(${PROJECT_NAME}-library-static PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNU}>)
#
#elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
#        ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") OR
#        ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
#    target_compile_options(${PROJECT_NAME}-library-shared PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNUCLANG}>)
#    target_compile_options(${PROJECT_NAME}-library-static PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNUCLANG}>)
#
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
#    target_compile_options(${PROJECT_NAME}-library-shared PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_MSVC}>)
#    target_compile_options(${PROJECT_NAME}-library-static PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_MSVC}>)
#    target_compile_options(${PROJECT_NAME}-library-shared PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS_MSVC}>)
#    target_compile_options(${PROJECT_NAME}-library-static PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS_MSVC}>)
#endif ()
#
##=======================================================================================================================
## Enable Link-Time Optimization
##=======================================================================================================================
#if (CMAKE_BUILD_TYPE STREQUAL "Release")
#    include(CheckIPOSupported)
#    check_ipo_supported(RESULT result OUTPUT output)
#    if (result)
#        set_property(TARGET ${PROJECT_NAME}-library-shared PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
#    endif ()
#endif ()



#=======================================================================================================================
# Install
#=======================================================================================================================
# Some basic stuff we'll need in this section
#include(GNUInstallDirs)
#include(CMakePackageConfigHelpers)
#set(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/openxlsx)
#
## Install interface headers
#install(
#        FILES ${OPENXLSX_CXX_INTERFACE_HEADERS}
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/openxlsx/${dir}
#)
#
## Install export header
#install(
#        FILES ${CMAKE_CURRENT_BINARY_DIR}/openxlsx_export.h
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/openxlsx
#)
#
## Targets
#install(
#        TARGETS OpenXLSX-static OpenXLSX-shared
#        EXPORT OpenXLSXTargets
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        COMPONENT lib
#        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        COMPONENT lib
#        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#        COMPONENT bin
#        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenXLSX
#)
#
## Package version
#write_basic_package_version_file(
#        "${CMAKE_CURRENT_BINARY_DIR}/OpenXLSX/OpenXLSXConfigVersion.cmake"
#        VERSION ${OpenXLSX_VERSION}
#        COMPATIBILITY AnyNewerVersion
#)
#install(
#        FILES
#        OpenXLSXConfig.cmake
#        "${CMAKE_CURRENT_BINARY_DIR}/OpenXLSX/OpenXLSXConfigVersion.cmake"
#        DESTINATION ${ConfigPackageLocation}
#)
#
## Package configuration
#configure_file(OpenXLSXConfig.cmake
#        "${CMAKE_CURRENT_BINARY_DIR}/OpenXLSX/OpenXLSXConfig.cmake"
#        COPYONLY
#        )
#
## Package export targets
#export(
#        EXPORT OpenXLSXTargets
#        FILE "${CMAKE_CURRENT_BINARY_DIR}/OpenXLSX/OpenXLSXTargets.cmake"
#        NAMESPACE OpenXLSX::
#)
#install(
#        EXPORT OpenXLSXTargets
#        FILE OpenXLSXTargets.cmake
#        NAMESPACE OpenXLSX::
#        DESTINATION ${ConfigPackageLocation}
#)



