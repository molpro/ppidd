cmake_minimum_required(VERSION 3.11)

project(ppidd LANGUAGES CXX)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(PPIDD_FORTRAN "Whether to build fortran binding" ON)
endif ()
if (PPIDD_FORTRAN)
    enable_language(Fortran)
endif ()


if (NOT MPI_FOUND)
    find_package(MPI REQUIRED)
endif (NOT MPI_FOUND)

if (CMAKE_Fortran_COMPILER_ID)
    message("Configuring Fortran support for ${CMAKE_CURRENT_SOURCE_DIR}")
    include(CheckFortranCompilerFlag)
    foreach (f "-fdefault-integer-8" "-i8")
        CHECK_Fortran_COMPILER_FLAG(${f} _fortran_flags)
        if (_fortran_flags)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${f}")
        endif ()
        unset(_fortran_flags CACHE)
    endforeach ()
endif ()


add_subdirectory(src)

enable_testing()
add_subdirectory(test)

