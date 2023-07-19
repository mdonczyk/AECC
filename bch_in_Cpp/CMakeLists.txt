cmake_minimum_required(VERSION 3.10)
project(BCH_code)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

set(BCH6351_SOURCE
     bch6351.cpp
     bch_logger.hpp)

# set(BCH6345_SOURCE
#     bch6345.cpp
#     bch_logger.hpp)

# set(BCH4836_SOURCE
#     bch4836.cpp
#     bch_logger.hpp)

# set(BCH4830_SOURCE
#     bch4830.cpp
#     bch_logger.hpp)

add_executable(bch6351 ${BCH6351_SOURCE})
# add_executable(bch6345 ${BCH6345_SOURCE})
# add_executable(bch4836 ${BCH4836_SOURCE})
# add_executable(bch4830 ${BCH4830_SOURCE})
