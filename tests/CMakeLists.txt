#=======================================================================================================================
# Download and define Catch2 library target
#=======================================================================================================================
add_library(Catch INTERFACE)
target_include_directories(Catch SYSTEM INTERFACE ${CMAKE_CURRENT_LIST_DIR}/Catch2)

#=======================================================================================================================
# Define TEST target
#=======================================================================================================================
add_executable(${PROJECT_NAME}-tests EXCLUDE_FROM_ALL "")
target_sources(${PROJECT_NAME}-tests
        PRIVATE
        TestMain.cpp
        PCComponent-test.cpp
        VPAntoineExt-test.cpp
        VPWagner-test.cpp
        VPRiedel-test.cpp
        VPAmbroseWalton-test.cpp
        VPHoffmannFlorin-test.cpp
        CDJoback-test.cpp
        SLVRackett-test.cpp
        SLVYenWoods-test.cpp
        SLVHankinsonThomson-test.cpp
        SLVElbro-test.cpp
        CLVThomson-test.cpp
        CLVAalto-test.cpp
        )

target_link_libraries(${PROJECT_NAME}-tests
        PRIVATE
        ${PROJECT_NAME}-library-static
        Catch
        )

#=======================================================================================================================
# Enable code coverage
#=======================================================================================================================
if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND NOT WIN32)
        target_compile_options(${PROJECT_NAME}-tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(${PROJECT_NAME}-tests PRIVATE -fprofile-instr-generate)
        #Uncomment in case of linker errors
        #link_libraries(clang_rt.profile-x86_64)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(${PROJECT_NAME}-tests PRIVATE --coverage)
        target_link_options(${PROJECT_NAME}-tests PRIVATE --coverage)
        #Uncomment in case of linker errors
        #link_libraries(gcov)
    endif ()
endif ()