cmake_minimum_required(VERSION 3.10)
project(BCH_code)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

set(BCH_SOURCE
     bch_simulator.cpp
     bch_utils.cpp
     bch_math.cpp
     bch_logger.hpp)

add_executable(bch_simulator ${BCH_SOURCE})
