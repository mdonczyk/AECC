# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mdonczyk/AECC/bch_in_Cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mdonczyk/AECC/bch_in_Cpp/build

# Include any dependencies generated for this target.
include CMakeFiles/bch4830xd.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/bch4830xd.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/bch4830xd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bch4830xd.dir/flags.make

CMakeFiles/bch4830xd.dir/bch4830.cpp.o: CMakeFiles/bch4830xd.dir/flags.make
CMakeFiles/bch4830xd.dir/bch4830.cpp.o: ../bch4830.cpp
CMakeFiles/bch4830xd.dir/bch4830.cpp.o: CMakeFiles/bch4830xd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mdonczyk/AECC/bch_in_Cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/bch4830xd.dir/bch4830.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/bch4830xd.dir/bch4830.cpp.o -MF CMakeFiles/bch4830xd.dir/bch4830.cpp.o.d -o CMakeFiles/bch4830xd.dir/bch4830.cpp.o -c /home/mdonczyk/AECC/bch_in_Cpp/bch4830.cpp

CMakeFiles/bch4830xd.dir/bch4830.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bch4830xd.dir/bch4830.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mdonczyk/AECC/bch_in_Cpp/bch4830.cpp > CMakeFiles/bch4830xd.dir/bch4830.cpp.i

CMakeFiles/bch4830xd.dir/bch4830.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bch4830xd.dir/bch4830.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mdonczyk/AECC/bch_in_Cpp/bch4830.cpp -o CMakeFiles/bch4830xd.dir/bch4830.cpp.s

# Object files for target bch4830xd
bch4830xd_OBJECTS = \
"CMakeFiles/bch4830xd.dir/bch4830.cpp.o"

# External object files for target bch4830xd
bch4830xd_EXTERNAL_OBJECTS =

../bch4830xd: CMakeFiles/bch4830xd.dir/bch4830.cpp.o
../bch4830xd: CMakeFiles/bch4830xd.dir/build.make
../bch4830xd: CMakeFiles/bch4830xd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mdonczyk/AECC/bch_in_Cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bch4830xd"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bch4830xd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/bch4830xd.dir/build: ../bch4830xd
.PHONY : CMakeFiles/bch4830xd.dir/build

CMakeFiles/bch4830xd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bch4830xd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bch4830xd.dir/clean

CMakeFiles/bch4830xd.dir/depend:
	cd /home/mdonczyk/AECC/bch_in_Cpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mdonczyk/AECC/bch_in_Cpp /home/mdonczyk/AECC/bch_in_Cpp /home/mdonczyk/AECC/bch_in_Cpp/build /home/mdonczyk/AECC/bch_in_Cpp/build /home/mdonczyk/AECC/bch_in_Cpp/build/CMakeFiles/bch4830xd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bch4830xd.dir/depend

