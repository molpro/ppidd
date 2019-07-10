cmake_minimum_required(VERSION 3.14)

project(ppidd LANGUAGES CXX)

option(PPIDD_FORTRAN "Whether to build Fortran binding" OFF)
option(INTEGER8 "Whether to compile with 8-byte (rather than 4-byte) Fortran integers" ON)
option(PPIDD-GA "Whether to build PPIDD on Global Arrays instead of pure MPI" ON)

if (PPIDD_FORTRAN)
    enable_language(Fortran)
endif ()

add_subdirectory(src)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    enable_testing()
    add_subdirectory(test)
endif ()
