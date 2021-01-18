#=======================================================================================================================
# Define TEST target
#=======================================================================================================================
add_executable(PClib-tests EXCLUDE_FROM_ALL "")
target_sources(PClib-tests
        PRIVATE
        TestMain.cpp
        Component/PCComponent-test.cpp
        Constants/CDJoback-test.cpp
        VaporPressure/VPAntoineExt-test.cpp
        VaporPressure/VPWagner-test.cpp
        VaporPressure/VPRiedel-test.cpp
        VaporPressure/VPAmbroseWalton-test.cpp
        VaporPressure/VPHoffmannFlorin-test.cpp
        LiquidVolume/SLVRackett-test.cpp
        LiquidVolume/SLVYenWoods-test.cpp
        LiquidVolume/SLVHankinsonThomson-test.cpp
        LiquidVolume/SLVElbro-test.cpp
        LiquidVolume/CLVThomson-test.cpp
        LiquidVolume/CLVAalto-test.cpp
        )

target_link_libraries(PClib-tests PRIVATE dependencies pclib)