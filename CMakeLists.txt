cmake_minimum_required(VERSION 3.17)
project(untitled)

set(CMAKE_CXX_STANDARD 14)

add_executable(untitled main.cpp Node.cpp Node.h Element.cpp Element.h GlobalData.cpp GlobalData.h Functions.cpp Functions.h ElemSolve.cpp ElemSolve.h)