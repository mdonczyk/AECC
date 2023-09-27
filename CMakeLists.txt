cmake_minimum_required(VERSION 3.10)
project(BCH_code)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

set(BCH_SOURCE
     main.cpp
     bch_simulator.cpp
     bch_utils.cpp
     bch_math.cpp
     bch_logger.hpp)

add_executable(bch_simulator ${BCH_SOURCE})

# Set warning level and warning flags
target_compile_options(bch_simulator PRIVATE
	-Wall                   # Enable all commonly used warnings
        -Wextra                 # Enable extra warnings
        -Wpedantic              # Enable strict ISO C++ compliance
        -Werror                 # Treat all warnings as errors
	-Wconversion            # Warn on implicit type conversions
	-Wsign-conversion       # Warn on sign conversion
        -Wformat=2              # Check format string arguments
        -Wshadow                # Warn when a local variable shadows another
        -Wundef                 # Warn if a preprocessor directive evaluates to undefined
        -Wnull-dereference      # Warn on potential null pointer dereference
        -Wdouble-promotion      # Warn on implicit float/double promotion
        -Wformat-signedness     # Warn on format string mismatches for signed/unsigned types
	-Wcast-align            # Warn on pointer casting that may result in incorrect alignment
)

