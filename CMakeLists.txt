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

# Eigen headers (vendored or fetched). Needed for optional 2Q synthesis.
option(QASMTRANS_USE_EIGEN "Enable Eigen-based features." OFF)
option(QASMTRANS_FETCH_EIGEN "Download Eigen if not vendored." OFF)
set(QASMTRANS_EIGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen-3.4.1"
    CACHE PATH "Path to Eigen headers (directory containing Eigen/)")

if(QASMTRANS_USE_EIGEN)
    if(EXISTS "${QASMTRANS_EIGEN_DIR}/Eigen/Core")
        set(EIGEN3_INCLUDE_DIR "${QASMTRANS_EIGEN_DIR}")
    else()
        if(QASMTRANS_FETCH_EIGEN)
            include(FetchContent)
            FetchContent_Declare(
                eigen3
                URL https://gitlab.com/libeigen/eigen/-/archive/3.4.1/eigen-3.4.1.tar.bz2
            )
            FetchContent_MakeAvailable(eigen3)
            set(EIGEN3_INCLUDE_DIR "${eigen3_SOURCE_DIR}")
        else()
            message(FATAL_ERROR "Eigen not found. Set QASMTRANS_EIGEN_DIR or enable QASMTRANS_FETCH_EIGEN=ON.")
        endif()
    endif()
endif()

# Add executable
add_executable(QASMTrans
    src/qasmtrans.cpp
    src/cli_support.cpp
    src/qick_emitter.cpp
)

# Make the project's headers and the vendored Lemon headers visible to the target
target_include_directories(QASMTrans
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party
)
if(QASMTRANS_USE_EIGEN)
    target_include_directories(QASMTrans PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_compile_definitions(QASMTrans PRIVATE QASMTRANS_USE_EIGEN=1)
endif()

# Link pybind
target_link_libraries(QASMTrans PRIVATE pybind11::embed Python3::Python)

# Python binding skeleton (pybind11 module)
pybind11_add_module(qasmtrans_core MODULE
    src/python_binding.cpp
    src/cli_support.cpp
)
target_include_directories(qasmtrans_core
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party
)
if(QASMTRANS_USE_EIGEN)
    target_include_directories(qasmtrans_core PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_compile_definitions(qasmtrans_core PRIVATE QASMTRANS_USE_EIGEN=1)
endif()
# Place the module next to the python/ helpers for easy import during development
set_target_properties(qasmtrans_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python
)

# Apple-provided Clang (and the Homebrew LLVM toolchain that scikit-build often
# selects when building wheels) does not always add the matching libc++/c++abi
# search paths when linking executables under pip. That leads to link failures
# such as "std::__1::__hash_memory" being undefined when QASMTrans is built
# during `pip install .`. Locate the runtime libraries next to the active
# compiler and link against them explicitly so we always pull in the matching
# C++ runtime regardless of the surrounding environment.
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    get_filename_component(_cxx_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(_toolchain_root "${_cxx_dir}/.." ABSOLUTE)
    set(_cxx_runtime_hint "${_toolchain_root}/lib/c++")
    set(_cxx_lib_hints "${_cxx_runtime_hint};${_toolchain_root}/lib")
    find_library(_libcxx c++ HINTS ${_cxx_lib_hints})
    find_library(_libcxxabi c++abi HINTS ${_cxx_lib_hints})
    foreach(target_name IN ITEMS QASMTrans qasmtrans_core)
        if(TARGET ${target_name})
            if(_libcxx)
                target_link_libraries(${target_name} PRIVATE ${_libcxx})
            endif()
            if(_libcxxabi)
                target_link_libraries(${target_name} PRIVATE ${_libcxxabi})
            endif()
        endif()
    endforeach()
endif()

install(TARGETS qasmtrans_core
    LIBRARY DESTINATION qasmtrans
    RUNTIME DESTINATION qasmtrans
    ARCHIVE DESTINATION qasmtrans
)
install(FILES
    python/qasmtrans/__init__.py
    python/qasmtrans/plot_pulses.py
    python/qasmtrans/qick_pulse_runner.py
    DESTINATION qasmtrans
)
