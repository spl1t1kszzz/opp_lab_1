cmake_minimum_required(VERSION 3.24)
project(opp_lab_1 C)

set(CMAKE_C_STANDARD 11)
include_directories(/opt/homebrew/bin/mpi)
add_executable(opp_lab_1
        main.c)
