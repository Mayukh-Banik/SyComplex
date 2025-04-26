cmake_minimum_required(VERSION 4.0)

project(syDTypes)

add_library(syDTypes INTERFACE)

target_include_directories(syDTypes INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/syDTypes
)
