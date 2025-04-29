set(EXECUTABLE_NAME complexTests)

add_executable(${EXECUTABLE_NAME} 
    tests.cpp
)

target_link_libraries(${EXECUTABLE_NAME}
    PRIVATE
        gtest
        gtest_main
        syComplex
)
