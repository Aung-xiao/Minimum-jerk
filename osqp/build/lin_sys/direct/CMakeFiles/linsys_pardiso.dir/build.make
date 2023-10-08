# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build

# Include any dependencies generated for this target.
include lin_sys/direct/CMakeFiles/linsys_pardiso.dir/depend.make

# Include the progress variables for this target.
include lin_sys/direct/CMakeFiles/linsys_pardiso.dir/progress.make

# Include the compile flags for this target's objects.
include lin_sys/direct/CMakeFiles/linsys_pardiso.dir/flags.make

lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.o: lin_sys/direct/CMakeFiles/linsys_pardiso.dir/flags.make
lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.o: ../lin_sys/direct/pardiso/pardiso_interface.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.o"
	cd /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.o   -c /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/lin_sys/direct/pardiso/pardiso_interface.c

lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.i"
	cd /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/lin_sys/direct/pardiso/pardiso_interface.c > CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.i

lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.s"
	cd /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/lin_sys/direct/pardiso/pardiso_interface.c -o CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.s

lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.o: lin_sys/direct/CMakeFiles/linsys_pardiso.dir/flags.make
lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.o: ../lin_sys/direct/pardiso/pardiso_loader.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.o"
	cd /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.o   -c /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/lin_sys/direct/pardiso/pardiso_loader.c

lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.i"
	cd /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/lin_sys/direct/pardiso/pardiso_loader.c > CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.i

lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.s"
	cd /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/lin_sys/direct/pardiso/pardiso_loader.c -o CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.s

linsys_pardiso: lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_interface.c.o
linsys_pardiso: lin_sys/direct/CMakeFiles/linsys_pardiso.dir/pardiso/pardiso_loader.c.o
linsys_pardiso: lin_sys/direct/CMakeFiles/linsys_pardiso.dir/build.make

.PHONY : linsys_pardiso

# Rule to build all files generated by this target.
lin_sys/direct/CMakeFiles/linsys_pardiso.dir/build: linsys_pardiso

.PHONY : lin_sys/direct/CMakeFiles/linsys_pardiso.dir/build

lin_sys/direct/CMakeFiles/linsys_pardiso.dir/clean:
	cd /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct && $(CMAKE_COMMAND) -P CMakeFiles/linsys_pardiso.dir/cmake_clean.cmake
.PHONY : lin_sys/direct/CMakeFiles/linsys_pardiso.dir/clean

lin_sys/direct/CMakeFiles/linsys_pardiso.dir/depend:
	cd /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/lin_sys/direct /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct /home/spy/hg_ws/src/MPC-Application-in-Trajectory-Tracking/mpc_mec/osqp/build/lin_sys/direct/CMakeFiles/linsys_pardiso.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lin_sys/direct/CMakeFiles/linsys_pardiso.dir/depend
