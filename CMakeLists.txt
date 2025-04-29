cmake_minimum_required (VERSION 3.28)

project(syComplex CXX)

add_library(syComplex INTERFACE)

target_include_directories(syComplex INTERFACE ${CMAKE_SOURCE_DIR}/syComplex)

option(ENABLE_TESTING "Enable testing" OFF)

if(ENABLE_TESTING)
    set(GOOGLETEST_DIR ${CMAKE_SOURCE_DIR}/googletest)
    add_subdirectory(${GOOGLETEST_DIR} ${CMAKE_BINARY_DIR}/googletest)
    add_subdirectory(tests)
endif()