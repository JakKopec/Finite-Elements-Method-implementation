# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\Programy\CLion\CLion 2020.2.4\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\Programy\CLion\CLion 2020.2.4\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\Studia\Semestr5\Mes\untitled

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\Studia\Semestr5\Mes\untitled\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/untitled.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/untitled.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/untitled.dir/flags.make

CMakeFiles/untitled.dir/main.cpp.obj: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/untitled.dir/main.cpp.obj"
	D:\Programy\MinGW\mingw32\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\untitled.dir\main.cpp.obj -c D:\Studia\Semestr5\Mes\untitled\main.cpp

CMakeFiles/untitled.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/main.cpp.i"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Studia\Semestr5\Mes\untitled\main.cpp > CMakeFiles\untitled.dir\main.cpp.i

CMakeFiles/untitled.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/main.cpp.s"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Studia\Semestr5\Mes\untitled\main.cpp -o CMakeFiles\untitled.dir\main.cpp.s

CMakeFiles/untitled.dir/Node.cpp.obj: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Node.cpp.obj: ../Node.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/untitled.dir/Node.cpp.obj"
	D:\Programy\MinGW\mingw32\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\untitled.dir\Node.cpp.obj -c D:\Studia\Semestr5\Mes\untitled\Node.cpp

CMakeFiles/untitled.dir/Node.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Node.cpp.i"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Studia\Semestr5\Mes\untitled\Node.cpp > CMakeFiles\untitled.dir\Node.cpp.i

CMakeFiles/untitled.dir/Node.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Node.cpp.s"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Studia\Semestr5\Mes\untitled\Node.cpp -o CMakeFiles\untitled.dir\Node.cpp.s

CMakeFiles/untitled.dir/Element.cpp.obj: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Element.cpp.obj: ../Element.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/untitled.dir/Element.cpp.obj"
	D:\Programy\MinGW\mingw32\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\untitled.dir\Element.cpp.obj -c D:\Studia\Semestr5\Mes\untitled\Element.cpp

CMakeFiles/untitled.dir/Element.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Element.cpp.i"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Studia\Semestr5\Mes\untitled\Element.cpp > CMakeFiles\untitled.dir\Element.cpp.i

CMakeFiles/untitled.dir/Element.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Element.cpp.s"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Studia\Semestr5\Mes\untitled\Element.cpp -o CMakeFiles\untitled.dir\Element.cpp.s

CMakeFiles/untitled.dir/GlobalData.cpp.obj: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/GlobalData.cpp.obj: ../GlobalData.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/untitled.dir/GlobalData.cpp.obj"
	D:\Programy\MinGW\mingw32\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\untitled.dir\GlobalData.cpp.obj -c D:\Studia\Semestr5\Mes\untitled\GlobalData.cpp

CMakeFiles/untitled.dir/GlobalData.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/GlobalData.cpp.i"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Studia\Semestr5\Mes\untitled\GlobalData.cpp > CMakeFiles\untitled.dir\GlobalData.cpp.i

CMakeFiles/untitled.dir/GlobalData.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/GlobalData.cpp.s"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Studia\Semestr5\Mes\untitled\GlobalData.cpp -o CMakeFiles\untitled.dir\GlobalData.cpp.s

CMakeFiles/untitled.dir/Functions.cpp.obj: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Functions.cpp.obj: ../Functions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/untitled.dir/Functions.cpp.obj"
	D:\Programy\MinGW\mingw32\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\untitled.dir\Functions.cpp.obj -c D:\Studia\Semestr5\Mes\untitled\Functions.cpp

CMakeFiles/untitled.dir/Functions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Functions.cpp.i"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Studia\Semestr5\Mes\untitled\Functions.cpp > CMakeFiles\untitled.dir\Functions.cpp.i

CMakeFiles/untitled.dir/Functions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Functions.cpp.s"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Studia\Semestr5\Mes\untitled\Functions.cpp -o CMakeFiles\untitled.dir\Functions.cpp.s

CMakeFiles/untitled.dir/Elem2Solve.cpp.obj: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Elem2Solve.cpp.obj: ../Elem2Solve.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/untitled.dir/Elem2Solve.cpp.obj"
	D:\Programy\MinGW\mingw32\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\untitled.dir\Elem2Solve.cpp.obj -c D:\Studia\Semestr5\Mes\untitled\Elem2Solve.cpp

CMakeFiles/untitled.dir/Elem2Solve.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Elem2Solve.cpp.i"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Studia\Semestr5\Mes\untitled\Elem2Solve.cpp > CMakeFiles\untitled.dir\Elem2Solve.cpp.i

CMakeFiles/untitled.dir/Elem2Solve.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Elem2Solve.cpp.s"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Studia\Semestr5\Mes\untitled\Elem2Solve.cpp -o CMakeFiles\untitled.dir\Elem2Solve.cpp.s

CMakeFiles/untitled.dir/matrix.cpp.obj: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/matrix.cpp.obj: ../matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/untitled.dir/matrix.cpp.obj"
	D:\Programy\MinGW\mingw32\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\untitled.dir\matrix.cpp.obj -c D:\Studia\Semestr5\Mes\untitled\matrix.cpp

CMakeFiles/untitled.dir/matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/matrix.cpp.i"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Studia\Semestr5\Mes\untitled\matrix.cpp > CMakeFiles\untitled.dir\matrix.cpp.i

CMakeFiles/untitled.dir/matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/matrix.cpp.s"
	D:\Programy\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Studia\Semestr5\Mes\untitled\matrix.cpp -o CMakeFiles\untitled.dir\matrix.cpp.s

# Object files for target untitled
untitled_OBJECTS = \
"CMakeFiles/untitled.dir/main.cpp.obj" \
"CMakeFiles/untitled.dir/Node.cpp.obj" \
"CMakeFiles/untitled.dir/Element.cpp.obj" \
"CMakeFiles/untitled.dir/GlobalData.cpp.obj" \
"CMakeFiles/untitled.dir/Functions.cpp.obj" \
"CMakeFiles/untitled.dir/Elem2Solve.cpp.obj" \
"CMakeFiles/untitled.dir/matrix.cpp.obj"

# External object files for target untitled
untitled_EXTERNAL_OBJECTS =

untitled.exe: CMakeFiles/untitled.dir/main.cpp.obj
untitled.exe: CMakeFiles/untitled.dir/Node.cpp.obj
untitled.exe: CMakeFiles/untitled.dir/Element.cpp.obj
untitled.exe: CMakeFiles/untitled.dir/GlobalData.cpp.obj
untitled.exe: CMakeFiles/untitled.dir/Functions.cpp.obj
untitled.exe: CMakeFiles/untitled.dir/Elem2Solve.cpp.obj
untitled.exe: CMakeFiles/untitled.dir/matrix.cpp.obj
untitled.exe: CMakeFiles/untitled.dir/build.make
untitled.exe: CMakeFiles/untitled.dir/linklibs.rsp
untitled.exe: CMakeFiles/untitled.dir/objects1.rsp
untitled.exe: CMakeFiles/untitled.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable untitled.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\untitled.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/untitled.dir/build: untitled.exe

.PHONY : CMakeFiles/untitled.dir/build

CMakeFiles/untitled.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\untitled.dir\cmake_clean.cmake
.PHONY : CMakeFiles/untitled.dir/clean

CMakeFiles/untitled.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\Studia\Semestr5\Mes\untitled D:\Studia\Semestr5\Mes\untitled D:\Studia\Semestr5\Mes\untitled\cmake-build-debug D:\Studia\Semestr5\Mes\untitled\cmake-build-debug D:\Studia\Semestr5\Mes\untitled\cmake-build-debug\CMakeFiles\untitled.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/untitled.dir/depend

