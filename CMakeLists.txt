# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.20)

# Setup cc/CC
if(DEFINED ENV{cc} AND DEFINED ENV{CC})
    message("Environment variables cc and CC are set")
    set(CMAKE_C_COMPILER "$ENV{cc}")
    set(CMAKE_CXX_COMPILER "$ENV{CC}")
else()
    message("Environment variables cc and/or CC are not set, using default settings")
endif()

# Project's name
project(QASMTrans)

# Set C++17 as the standard for the entire project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the -O3 optimization flag for all configurations
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Set pybind environment stuff

set(Python3_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python" CACHE FILEPATH "")
set(Python3_ROOT_DIR "$ENV{CONDA_PREFIX}" CACHE PATH "")
set(pybind11_DIR "$ENV{CONDA_PREFIX}/share/cmake/pybind11" CACHE PATH "")
set(CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}" CACHE PATH "")

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Add executable
add_executable(QASMTrans
    src/qasmtrans.cpp
    src/qick_emitter.cpp
)


# Link pybind
target_link_libraries(QASMTrans PRIVATE pybind11::embed Python3::Python)
