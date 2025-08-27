cmake_minimum_required(VERSION 3.24.0-rc4)

set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(la_arbre_a_merde)

include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        v2.12.0
)

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

pybind11_add_module(${PROJECT_NAME} src/main.cpp)

add_subdirectory(src/hclust-cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE fastcluster)