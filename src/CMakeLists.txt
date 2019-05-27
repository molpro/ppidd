cmake_minimum_required(VERSION 3.14)
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
            #        GIT_REPOSITORY https://github.com/GlobalArrays/ga
            #        GIT_TAG v5.7
            GIT_REPOSITORY https://github.com/pjknowles/ga
            GIT_TAG c270f091
    )
    FetchContent_MakeAvailable(ga)
endif (ga-config)


set(PPIDD_DEFINES_H "#include \"ppidd_defines.h\"")
set(PPIDD_PROTOTYPES_H "#include \"ppidd_prototypes.h\"")
configure_file(ppidd.h.in ${CMAKE_CURRENT_BINARY_DIR}/ppidd.h)
add_library(ppidd
        "${CMAKE_CURRENT_BINARY_DIR}/ppidd.h" ppidd.cpp ppidd_defines.h
        ppidd_ga_mpi.h ppidd_ga_mpi.cpp
        ppidd_mpi2.h ppidd_mpi2.cpp
        ppidd_no_mpi.cpp
        ppidd_prototypes.h mpi_utils.cpp mpiga_base.cpp
        mpi_nxtval.cpp mpimutex-hybrid.cpp mpi_helpmutex.cpp)
if (CMAKE_Fortran_COMPILER_ID)
    configure_file(ppidd_module.F.in ${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F)
    target_sources(ppidd PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/ppidd_module.F")
endif ()

target_compile_features(ppidd PRIVATE cxx_constexpr)
target_compile_options(ppidd PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:AppleClang>>:-Wall>)
target_compile_definitions(ppidd PUBLIC FORTINT=int64_t HAVE_MPI_H HAVE_GA_H)
target_link_libraries(ppidd PUBLIC ga::ga MPI::MPI_CXX)
target_include_directories(ppidd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(ppidd PUBLIC ${CMAKE_CURRENT_BINARY_DIR}) # For Fortran modules
add_library(ppidd::ppidd ALIAS ppidd)

find_package(Doxygen)
if (DOXYGEN_FOUND)
    add_custom_target(${PROJECT_NAME}-doc ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            )
    get_property(SOURCES TARGET ppidd PROPERTY SOURCES)
    string(REPLACE ";/" " /" SOURCE1 "${SOURCES}")
    string(REPLACE ";" " ${CMAKE_CURRENT_SOURCE_DIR}/" SOURCE "${SOURCE1}")
    configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            DEPENDS ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif (DOXYGEN_FOUND)

