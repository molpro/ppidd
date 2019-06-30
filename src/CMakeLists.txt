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
        add_library(ga::ga INTERFACE IMPORTED)
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
if (PPIDD-GA)
    target_sources(ppidd PRIVATE ppidd_ga_mpi.h ppidd_ga_mpi.cpp)
    target_compile_definitions(ppidd PRIVATE HAVE_GA_H)
    target_link_libraries(ppidd PUBLIC ga::ga)
endif ()
if (CMAKE_Fortran_COMPILER_ID)
    configure_file(ppidd_module.F.in ${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F)
    target_sources(ppidd PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F")
endif ()

if (INTEGER8)
    target_compile_definitions(ppidd PUBLIC _I8_)
    target_compile_options(ppidd PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-fdefault-integer-8>)
    target_compile_options(ppidd PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:Intel>>:-i8>)
    target_compile_options(ppidd PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:PathScale>>:-i8>)
    target_compile_options(ppidd PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:PGI>>:-i8>)
    target_compile_options(ppidd PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:Cray>>:-s integer64>)
    target_compile_options(ppidd PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:SunPro>>:-xtypemap=integer:64>)
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
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
install(TARGETS ppidd EXPORT ppiddTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
        PUBLIC_HEADER DESTINATION include
        )
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

find_package(Doxygen)
if (DOXYGEN_FOUND)
    add_custom_target(${PROJECT_NAME}-doc ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            )
    get_property(SOURCES TARGET ppidd PROPERTY SOURCES)
    string(REPLACE ";/" " /" SOURCE1 "${SOURCES}")
    string(REPLACE ";" " ${CMAKE_CURRENT_SOURCE_DIR}/" SOURCE "${SOURCE1}")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            DEPENDS ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif (DOXYGEN_FOUND)

