#=======================================================================================================================
# EXTERNAL LIBRARIES
#   Define external libraries used by OpenXLSX. The libraries (Zippy and PugiXML) are header-only, so
#   INTERFACE libraries should be defined.
#=======================================================================================================================
add_library(external INTERFACE IMPORTED)
target_include_directories(external SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/external/>)


#add_library(Zippy INTERFACE IMPORTED)
#target_include_directories(Zippy SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/external/zippy/>)
#
#add_library(PugiXML INTERFACE IMPORTED)
#target_include_directories(PugiXML SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/external/pugixml/>)
#if (${ENABLE_COMPACT_MODE})
#    target_compile_definitions(PugiXML INTERFACE PUGIXML_COMPACT)
#endif ()
#
#add_library(NoWide INTERFACE IMPORTED)
#target_include_directories(NoWide SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/external/nowide/>)
#
#if (${ENABLE_UNICODE_FILENAMES})
#    add_compile_definitions(UNICODE_FILENAMES_ENABLED)
#endif ()


#=======================================================================================================================
# COMPILER FEATURES
#   Some older C++17 compilers don't support the char_conv features. If the compiler doesn't support it,
#   a less optimal work-around will be used.
#=======================================================================================================================
#include(CheckCXXSourceCompiles)
#check_cxx_source_compiles("
#                          #include <array>
#                          #include <string>
#                          #include <charconv>
#
#                          int main() {
#                                  std::array<char, 7> str {};
#                                  auto p = std::to_chars(str.data(), str.data() + str.size(), 12345).ptr;
#                                  auto strResult = std::string(str.data(), p - str.data());
#
#                                  unsigned long value = 0;
#                                  std::from_chars(strResult.data(), strResult.data() + strResult.size(), value);
#
#                                  return 0;
#                          }" CHARCONV_RESULT)
#
#if (CHARCONV_RESULT)
#    add_compile_definitions(CHARCONV_ENABLED)
#endif ()

#=======================================================================================================================
# PROJECT FILES
#   List of project source files
#=======================================================================================================================
set(PROJECT_TEMPLATE_SOURCES
        library/PCComponent.cpp
        library/VaporPressure/VPRiedel.cpp
        library/VaporPressure/VPAntoineExt.cpp
        library/VaporPressure/VPWagner.cpp
        library/VaporPressure/VPAmbroseWalton.cpp
        library/VaporPressure/VPHoffmannFlorin.cpp
        library/ConstantData/CDJoback.cpp
        library/LiquidVolume/SLVRackett.cpp
        library/LiquidVolume/SLVHankinsonThomson.cpp
        library/LiquidVolume/SLVYenWoods.cpp
        library/LiquidVolume/SLVElbro.cpp
        library/LiquidVolume/CLVThomson.cpp
        library/LiquidVolume/CLVAalto.cpp
        library/HeatCapacity/IGAlyLee.cpp
        library/EquationOfState/EOSPengRobinson.cpp
        )


#=======================================================================================================================
# STATIC LIBRARY
#   Define the static library
#=======================================================================================================================
add_library(${PROJECT_NAME}-library-static STATIC "")
target_sources(${PROJECT_NAME}-library-static PRIVATE ${PROJECT_TEMPLATE_SOURCES})
target_include_directories(${PROJECT_NAME}-library-static
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/library>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)     # For export header
target_link_libraries(${PROJECT_NAME}-library-static PRIVATE external)
target_compile_definitions(${PROJECT_NAME}-library-static PRIVATE OPENXLSX_STATIC_DEFINE)

#=======================================================================================================================
# SHARED LIBRARY
#   Define the shared library
#=======================================================================================================================
add_library(${PROJECT_NAME}-library-shared SHARED "")
target_sources(${PROJECT_NAME}-library-shared PRIVATE ${PROJECT_TEMPLATE_SOURCES})
target_include_directories(${PROJECT_NAME}-library-shared
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/library>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)     # For export header
target_link_libraries(${PROJECT_NAME}-library-shared PRIVATE external)

#=======================================================================================================================
# GENERATE EXPORT HEADER
#   Generate preprocessor definitions for setting symbol visibility.
#=======================================================================================================================
include(GenerateExportHeader)
generate_export_header(${PROJECT_NAME}-library-shared)
        #BASE_NAME openxlsx
        #EXPORT_FILE_NAME OpenXLSX-Exports.hpp
        #EXPORT_MACRO_NAME OPENXLSX_EXPORT
        #NO_EXPORT_MACRO_NAME OPENXLSX_HIDDEN)

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

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(${PROJECT_NAME}-library-shared PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNU}>)
    target_compile_options(${PROJECT_NAME}-library-static PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNU}>)

elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
        ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") OR
        ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    target_compile_options(${PROJECT_NAME}-library-shared PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNUCLANG}>)
    target_compile_options(${PROJECT_NAME}-library-static PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_GNUCLANG}>)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    target_compile_options(${PROJECT_NAME}-library-shared PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_MSVC}>)
    target_compile_options(${PROJECT_NAME}-library-static PRIVATE $<$<CONFIG:Debug>:${DEBUG_FLAGS_MSVC}>)
    target_compile_options(${PROJECT_NAME}-library-shared PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS_MSVC}>)
    target_compile_options(${PROJECT_NAME}-library-static PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS_MSVC}>)
endif ()

#=======================================================================================================================
# Enable Link-Time Optimization
#=======================================================================================================================
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output)
    if (result)
        set_property(TARGET ${PROJECT_NAME}-library-shared PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif ()
endif ()


#=======================================================================================================================
# Enable code coverage
#=======================================================================================================================
if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND NOT WIN32)
        target_compile_options(${PROJECT_NAME}-library-shared PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(${PROJECT_NAME}-library-shared PRIVATE -fprofile-instr-generate)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(${PROJECT_NAME}-library-shared PRIVATE --coverage)
        target_link_options(${PROJECT_NAME}-library-shared PRIVATE --coverage)
    endif ()
endif ()


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



