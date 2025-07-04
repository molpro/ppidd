cmake_minimum_required(VERSION 3.14)
if (NOT MPI_FOUND)
    find_package(MPI REQUIRED)
endif (NOT MPI_FOUND)
if (PPIDD-GA)
    message(STATUS "Building PPIDD on Global Arrays")
    set(PPIDD_IMPL_DEFAULT 2)
    if (NOT TARGET ga::ga)
        find_program(ga-config ga-config)
        if (ga-config AND NOT GA_BUILD)
            message(STATUS "Using installed Global Arrays")
            foreach (f cflags cppflags ldflags libs)
                execute_process(COMMAND ${ga-config} --${f} OUTPUT_VARIABLE ga-${f})
                string(STRIP "${ga-${f}}" ga-${f})
            endforeach ()
            add_library(ga::ga INTERFACE IMPORTED GLOBAL)
            target_link_libraries(ga::ga INTERFACE ${ga-libs})
            target_link_options(ga::ga INTERFACE ${ga-ldflags})
            target_compile_options(ga::ga INTERFACE ${ga-cflags} ${ga-cppflags})
        else (ga-config)
            message(STATUS "Building Global Arrays with GA_RUNTIME=${GA_RUNTIME}.")
            include(FetchContent)
            FetchContent_Declare(
                    ga
                    GIT_REPOSITORY https://github.com/GlobalArrays/ga
                    #                                        GIT_TAG 2518e23433385bfa3726d507b8cd0d7ed038021b
                    GIT_TAG v5.8.1
                    #                    GIT_TAG master
            )
            set(ENABLE_FORTRAN OFF)
            set(ENABLE_BLAS OFF)
            FetchContent_MakeAvailable(ga)
        endif ()
    endif ()
else ()
    message("Building PPIDD on pure MPI")
    set(PPIDD_IMPL_DEFAULT 3)
endif ()


set(PPIDD_DEFINES_H "#include \"ppidd_defines.h\"")
set(PPIDD_PROTOTYPES_H "#include \"ppidd_prototypes.h\"")
configure_file(ppidd.h.in ${CMAKE_CURRENT_BINARY_DIR}/ppidd.h)
LibraryManager_Add(ppidd NAMESPACE ppidd)
LibraryManager_Append(ppidd
        SOURCES
        ppidd.cpp
        ppidd_mpi2.h ppidd_mpi2.cpp
        ppidd_no_mpi.cpp
        mpi_utils.cpp mpiga_base.cpp
        mpi_nxtval.cpp mpimutex-hybrid.cpp mpi_helpmutex.cpp
        PUBLIC_HEADER
        ppidd_defines.h ppidd_prototypes.h mpiga_base.h mpimutex.h
        )
LibraryManager_AppendExternal(ppidd
        PUBLIC_HEADER
        "${CMAKE_CURRENT_BINARY_DIR}/ppidd.h"
        INCLUDE_DIR
        "${CMAKE_CURRENT_BINARY_DIR}"
        )
target_compile_definitions(ppidd PRIVATE BCAST_BATCH_SIZE=134217728)
if (PPIDD-GA)
    LibraryManager_Append(ppidd SOURCES ppidd_ga_mpi.h ppidd_ga_mpi.cpp)
    target_compile_definitions(ppidd PRIVATE HAVE_GA_H)
    if (NOT TARGET ga::ga)
        add_library(ga::ga ALIAS ga)
        include(LibraryManagerUtils)
        print_target_properties(ga)
        get_target_property(_ga-src ga SOURCE_DIR)
        get_target_property(_ga-binary ga BINARY_DIR)
        message(STATUS "ga::ga target not found, ga SOURCE_DIR=${_ga-src} BINARY_DIR=${_ga-binary}")
        target_include_directories(ga PUBLIC "$<BUILD_INTERFACE:${_ga-src}/global/src>")
        target_include_directories(ga PUBLIC "$<BUILD_INTERFACE:${_ga-src}/ma>")
        target_include_directories(ga PUBLIC "$<BUILD_INTERFACE:${_ga-binary}/ma>")
        target_include_directories(ga PUBLIC "$<BUILD_INTERFACE:${_ga-binary}/gaf2c>")
        target_include_directories(ga PUBLIC "$<BUILD_INTERFACE:${_ga-src}/pario/eaf>")
        target_include_directories(ga PUBLIC "$<BUILD_INTERFACE:${_ga-src}/pario/sf>")
    else ()
        message(STATUS "ga::ga target found")
    endif ()
    target_link_libraries(ppidd PUBLIC ga::ga)
endif ()
if (FORTRAN)
    configure_file(ppidd_module.F90.in ${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F90)
    message(DEBUG "Adding ppidd_module.F90")
    #    LibraryManager_AppendExternal(ppidd
    #            SOURCES "${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F90"
    #            INCLUDE_DIR
    #            "${CMAKE_CURRENT_BINARY_DIR}"
    #            )
    LibraryManager_Append(ppidd SOURCES
            SOURCES "${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F90"
            )
else ()
    message(DEBUG "Not building Fortran")
endif ()

target_compile_features(ppidd PRIVATE cxx_std_17)
target_compile_features(ppidd PRIVATE cxx_constexpr)
target_compile_options(ppidd PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:AppleClang>>:-Wall>)
target_compile_definitions(ppidd PRIVATE HAVE_MPI_H)
target_link_libraries(ppidd PUBLIC MPI::MPI_CXX)
#target_include_directories(ppidd PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
#        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
#        $<INSTALL_INTERFACE:include>
#        )
#add_library(ppidd::ppidd ALIAS ppidd)
#if (CMAKE_Fortran_MODULE_DIRECTORY)
#    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
#endif ()
#install(TARGETS ppidd EXPORT ppiddTargets
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib
#        INCLUDES DESTINATION include
#        PUBLIC_HEADER DESTINATION include
#        )
#install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ppidd.h" ppidd_defines.h ppidd_prototypes.h
#        DESTINATION include
#        )
if (PPIDD-GA)
    if (NOT ga-config)
        foreach (t ga ga_src dra eaf elio sf linalg)
            install(TARGETS ${t} EXPORT gaTargets
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
                    INCLUDES DESTINATION include
                    PUBLIC_HEADER DESTINATION include
                    )
        endforeach ()
        install(EXPORT gaTargets
                FILE gaTargets.cmake
                NAMESPACE ga::
                DESTINATION lib/cmake/ga)
    endif ()
endif ()
#install(EXPORT ppidd-export
#        FILE ppiddTargets.cmake
#        NAMESPACE ppidd::
#        DESTINATION lib/cmake/ppidd)
include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ppiddConfig.cmake"
        "include CMakeFindDependencyMacro)
        find_dependency(MPI)
        include(\"\${CMAKE_CURRENT_LIST_DIR}/ppiddTargets.cmake\")
        ")
write_basic_package_version_file("ppiddConfigVersion.cmake" VERSION 1.0 COMPATIBILITY SameMajorVersion)
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/ppiddConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/ppiddConfigVersion.cmake" DESTINATION lib/cmake/ppidd)

LibraryManager_Install(ppidd PKG_CONFIG)
LibraryManager_Export(ppidd)
