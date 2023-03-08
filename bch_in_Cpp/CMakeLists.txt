cmake_minimum_required(VERSION 3.10)
project(BCH_code)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

add_executable(bch6351 bch6351.cpp bch6351.hpp)

add_executable(bch6345 bch6345.cpp bch6345.hpp)

add_executable(bch4836 bch4836.cpp bch4836.hpp)

add_executable(bch4830 bch4830.cpp bch4830.hpp)
