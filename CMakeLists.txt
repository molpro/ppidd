cmake_minimum_required(VERSION 3.15)

project(ppidd LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG master
)
FetchContent_MakeAvailable(LibraryManager)

LibraryManager_Project(FORTRAN_OPTION MPI_OPTION)

option(PPIDD-GA "Whether to build PPIDD on Global Arrays instead of pure MPI" ON)
set(GA_RUNTIME "MPI_2SIDED" CACHE STRING "Global Arrays MPI runtime. Must be MPI_2SIDED, MPI_PROGRESS_RANK, MPI_RMA, MPI_MULTITHREADED, or MPI_PROGRESS_THREAD")
option(GA_BUILD "Whether to explicitly build Global Arrays instead of accepting system-installed library" OFF)

add_subdirectory(src)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    enable_testing()
    add_subdirectory(test)
endif ()
