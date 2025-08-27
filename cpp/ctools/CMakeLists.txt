cmake_minimum_required(VERSION 3.24.0-rc4)
project(ctools)

include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        v2.13.6
)

# FetchContent_MakeAvailable handles both downloading and configuring the dependency.
# It replaces both FetchContent_Populate() and add_subdirectory().
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(${PROJECT_NAME} src/main.cpp)