cmake_minimum_required(VERSION 3.10)
project(ME553_2024)

set(CMAKE_CXX_STANDARD 17)
include(cmake/helper.cmake)

# exercise

# examples
create_executable(rotation_example_angle_axis src/rotation_example_angle_axis.cpp)
create_executable(rotation_example_euler_angles src/rotation_example_euler_angles.cpp)
create_executable(panda_example src/panda.cpp)

# exam



