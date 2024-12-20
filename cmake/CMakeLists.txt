cmake_minimum_required(VERSION 3.22)
project(untitled)

set(CMAKE_CXX_STANDARD 20)

set(TYPES "FIXED(32,16), DOUBLE" CACHE STRING "Types to precompile")

add_executable(untitled main.cpp)
