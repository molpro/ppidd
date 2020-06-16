cmake_minimum_required(VERSION 3.14)
if (NOT MPI_FOUND)
    find_package(MPI REQUIRED)
endif (NOT MPI_FOUND)
if (PPIDD-GA)
    message("Building PPIDD on Global Arrays")
    set(PPIDD_IMPL_DEFAULT 2)
    find_program(ga-config ga-config)
    if (ga-config)
        foreach (f cflags cppflags ldflags libs)
            execute_process(COMMAND ${ga-config} --${f} OUTPUT_VARIABLE ga-${f})
            string(STRIP "${ga-${f}}" ga-${f})
        endforeach ()
        add_library(ga::ga INTERFACE IMPORTED GLOBAL)
        target_link_libraries(ga::ga INTERFACE ${ga-libs})
        target_link_options(ga::ga INTERFACE ${ga-ldflags})
        target_compile_options(ga::ga INTERFACE ${ga-cflags} ${ga-cppflags})
    else (ga-config)
        include(FetchContent)
        FetchContent_Declare(
                ga
                GIT_REPOSITORY https://github.com/GlobalArrays/ga
                GIT_TAG 2518e23433385bfa3726d507b8cd0d7ed038021b
                #                GIT_REPOSITORY https://github.com/pjknowles/ga
                #                GIT_TAG c270f091
        )
        if (PPIDD_FORTRAN)
            set(ENABLE_FORTRAN ON)
        endif ()
        FetchContent_MakeAvailable(ga)
    endif (ga-config)
else ()
    message("Building PPIDD on pure MPI")
    set(PPIDD_IMPL_DEFAULT 3)
endif ()


set(PPIDD_DEFINES_H "#include \"ppidd_defines.h\"")
set(PPIDD_PROTOTYPES_H "#include \"ppidd_prototypes.h\"")
configure_file(ppidd.h.in ${CMAKE_CURRENT_BINARY_DIR}/ppidd.h)
add_library(ppidd
        "${CMAKE_CURRENT_BINARY_DIR}/ppidd.h" ppidd.cpp ppidd_defines.h
        ppidd_mpi2.h ppidd_mpi2.cpp
        ppidd_no_mpi.cpp
        ppidd_prototypes.h mpi_utils.cpp mpiga_base.cpp
        mpi_nxtval.cpp mpimutex-hybrid.cpp mpi_helpmutex.cpp)
target_compile_definitions(ppidd PRIVATE BCAST_BATCH_SIZE=134217728)
if (PPIDD-GA)
    target_sources(ppidd PRIVATE ppidd_ga_mpi.h ppidd_ga_mpi.cpp)
    target_compile_definitions(ppidd PRIVATE HAVE_GA_H)
    target_link_libraries(ppidd PUBLIC ga::ga)
endif ()
if (CMAKE_Fortran_COMPILER_ID)
    configure_file(ppidd_module.F90.in ${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F90)
    target_sources(ppidd PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F90")
endif ()

if (INTEGER8)
    target_compile_definitions(ppidd PRIVATE $<$<COMPILE_LANGUAGE:CXX>:FORTINT=int64_t>)
else ()
    target_compile_definitions(ppidd PRIVATE $<$<COMPILE_LANGUAGE:CXX>:FORTINT=int32_t>)
endif ()

target_compile_features(ppidd PRIVATE cxx_constexpr)
target_compile_options(ppidd PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:AppleClang>>:-Wall>)
target_compile_definitions(ppidd PRIVATE HAVE_MPI_H)
target_link_libraries(ppidd PUBLIC MPI::MPI_CXX)
target_include_directories(ppidd PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
        $<INSTALL_INTERFACE:include>
        )
add_library(ppidd::ppidd ALIAS ppidd)
if (CMAKE_Fortran_MODULE_DIRECTORY)
    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
endif ()
install(TARGETS ppidd EXPORT ppiddTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
        PUBLIC_HEADER DESTINATION include
        )
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ppidd.h" ppidd_defines.h ppidd_prototypes.h
        DESTINATION include
        )
if (PPIDD-GA)
    if (NOT ga-config)
        foreach (t ga ga_src dra eaf elio sf linalg)
            message("Export target ${t}")
            install(TARGETS ${t} EXPORT ppiddTargets
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
                    INCLUDES DESTINATION include
                    PUBLIC_HEADER DESTINATION include
                    )
        endforeach ()
    endif ()
endif ()
install(EXPORT ppiddTargets
        FILE ppiddTargets.cmake
        NAMESPACE ppidd::
        DESTINATION lib/cmake/ppidd)
include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ppiddConfig.cmake"
        "include CMakeFindDependencyMacro)
        find_dependency(MPI)
        include(\"\${CMAKE_CURRENT_LIST_DIR}/ppiddTargets.cmake\")
        ")
write_basic_package_version_file("ppiddConfigVersion.cmake" VERSION 1.0 COMPATIBILITY SameMajorVersion)
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/ppiddConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/ppiddConfigVersion.cmake" DESTINATION lib/cmake/ppidd)
