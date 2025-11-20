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

set(Python3_FIND_VIRTUALENV FIRST)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE _pybind11_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pybind11_status
)

if (NOT _pybind11_status EQUAL 0)
    message(FATAL_ERROR "Failed to find pybind11 via Python. Did you run 'pip install pybind11'?")
endif()

list(APPEND CMAKE_PREFIX_PATH "${_pybind11_dir}")

find_package(pybind11 CONFIG REQUIRED)

# Add executable
add_executable(QASMTrans
    src/qasmtrans.cpp
    src/qick_emitter.cpp
)

# Make the project's headers and the vendored Lemon headers visible to the target
target_include_directories(QASMTrans
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party
)

# Link pybind
target_link_libraries(QASMTrans PRIVATE pybind11::embed Python3::Python)
