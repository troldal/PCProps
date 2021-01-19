#=======================================================================================================================
# Define TEST target
#=======================================================================================================================
add_executable(EOSlib-tests EXCLUDE_FROM_ALL "")
target_sources(EOSlib-tests
        PRIVATE
        TestMain.cpp
        EOSPengRobinson-test.cpp
        )

target_link_libraries(EOSlib-tests PRIVATE dependencies PropertyLib EOSLib)