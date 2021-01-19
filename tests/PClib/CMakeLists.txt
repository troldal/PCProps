#=======================================================================================================================
# Define TEST target
#=======================================================================================================================
add_executable(PClib-tests EXCLUDE_FROM_ALL "")
target_sources(PClib-tests
        PRIVATE
        TestMain.cpp
        Component/PCComponent-test.cpp
        Constants/Joback-test.cpp
        VaporPressure/AntoineExt-test.cpp
        VaporPressure/Wagner-test.cpp
        VaporPressure/Riedel-test.cpp
        VaporPressure/AmbroseWalton-test.cpp
        VaporPressure/HoffmannFlorin-test.cpp
        SaturatedLiquidVolume/Rackett-test.cpp
        SaturatedLiquidVolume/YenWoods-test.cpp
        SaturatedLiquidVolume/HankinsonThomson-test.cpp
        SaturatedLiquidVolume/Elbro-test.cpp
        CompressedLiquidVolume/Thomson-test.cpp
        CompressedLiquidVolume/Aalto-test.cpp
        )

target_link_libraries(PClib-tests PRIVATE dependencies PropertyLib)