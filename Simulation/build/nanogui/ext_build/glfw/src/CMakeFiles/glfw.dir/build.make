# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vismay/ElasticBodies/Simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vismay/ElasticBodies/Simulation/build

# Include any dependencies generated for this target.
include nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/depend.make

# Include the progress variables for this target.
include nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/progress.make

# Include the compile flags for this target's objects.
include nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/flags.make

# Object files for target glfw
glfw_OBJECTS =

# External object files for target glfw
glfw_EXTERNAL_OBJECTS = \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/context.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/init.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/input.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/monitor.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/window.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/x11_init.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/x11_monitor.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/x11_window.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/xkb_unicode.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/linux_joystick.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_time.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_tls.c.o" \
"/home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/glx_context.c.o"

nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/context.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/init.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/input.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/monitor.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/window.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/x11_init.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/x11_monitor.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/x11_window.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/xkb_unicode.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/linux_joystick.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_time.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_tls.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/glx_context.c.o
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/build.make
nanogui/ext_build/glfw/src/libglfw3.a: nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vismay/ElasticBodies/Simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking C static library libglfw3.a"
	cd /home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src && $(CMAKE_COMMAND) -P CMakeFiles/glfw.dir/cmake_clean_target.cmake
	cd /home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/glfw.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/build: nanogui/ext_build/glfw/src/libglfw3.a

.PHONY : nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/build

nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/requires:

.PHONY : nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/requires

nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/clean:
	cd /home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src && $(CMAKE_COMMAND) -P CMakeFiles/glfw.dir/cmake_clean.cmake
.PHONY : nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/clean

nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/depend:
	cd /home/vismay/ElasticBodies/Simulation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vismay/ElasticBodies/Simulation /home/vismay/libigl/external/nanogui/ext/glfw/src /home/vismay/ElasticBodies/Simulation/build /home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src /home/vismay/ElasticBodies/Simulation/build/nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : nanogui/ext_build/glfw/src/CMakeFiles/glfw.dir/depend

