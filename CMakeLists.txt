cmake_minimum_required(VERSION 3.10)
project(ME553_2025)

set(CMAKE_CXX_STANDARD 17)
include(cmake/helper.cmake)

add_subdirectory(third_party/raisimLib)

# exercise
create_executable(exercise1 src/exercise1.cpp)
create_executable(exercise2 src/exercise2.cpp)

# examples
create_executable(rotation_example_angle_axis src/rotation_example_angle_axis.cpp)
create_executable(rotation_example_euler_angles src/rotation_example_euler_angles.cpp)
create_executable(panda_example src/panda.cpp)

# exam
