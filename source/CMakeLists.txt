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


#=======================================================================================================================
# LOAD SUBDIRECTORIES
#   ---
#=======================================================================================================================

add_subdirectory(common)
add_subdirectory(Interfaces)
add_subdirectory(Properties)
add_subdirectory(EquationOfState)
add_subdirectory(Correlations)
add_subdirectory(PureComponent)
add_subdirectory(Fluid)
add_subdirectory(UnitOps)
add_subdirectory(Numerics)


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



