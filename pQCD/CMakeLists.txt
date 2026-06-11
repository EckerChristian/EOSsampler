cmake_minimum_required(VERSION 3.10)
project(pQCD)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(PQCD pQCD.cpp definitions.h)