cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

include_directories(BEFORE pylibseq/src/libsequence)
add_subdirectory(pylibseq)

install(TARGETS _libsequence DESTINATION ${SKBUILD_PROJECT_NAME})
