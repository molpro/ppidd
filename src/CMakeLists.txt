cmake_minimum_required(VERSION 3.11)
set(PPIDD_IMPL_DEFAULT 2)

find_program(ga-config ga-config)
if (ga-config)
    foreach (f cflags cppflags ldflags libs)
        execute_process(COMMAND ${ga-config} --${f} OUTPUT_VARIABLE ga-${f})
        string(STRIP "${ga-${f}}" ga-${f})
        message(ga-${f}:${ga-${f}})
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
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
            DEPENDS ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif (DOXYGEN_FOUND)

# debugging below here
# Get all propreties that cmake supports
execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)

# Convert command output into a CMake list
STRING(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
STRING(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
# Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
list(FILTER CMAKE_PROPERTY_LIST EXCLUDE REGEX "^LOCATION$|^LOCATION_|_LOCATION$")
# For some reason, "TYPE" shows up twice - others might too?
list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)

# build whitelist by filtering down from CMAKE_PROPERTY_LIST in case cmake is
# a different version, and one of our hardcoded whitelisted properties
# doesn't exist!
unset(CMAKE_WHITELISTED_PROPERTY_LIST)
foreach(prop ${CMAKE_PROPERTY_LIST})
    if(prop MATCHES "^(INTERFACE|[_a-z]|IMPORTED_LIBNAME_|MAP_IMPORTED_CONFIG_)|^(COMPATIBLE_INTERFACE_(BOOL|NUMBER_MAX|NUMBER_MIN|STRING)|EXPORT_NAME|IMPORTED(_GLOBAL|_CONFIGURATIONS|_LIBNAME)?|NAME|TYPE|NO_SYSTEM_FROM_IMPORTED)$")
        list(APPEND CMAKE_WHITELISTED_PROPERTY_LIST ${prop})
    endif()
endforeach(prop)

function(print_properties)
    message ("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
endfunction(print_properties)

function(print_whitelisted_properties)
    message ("CMAKE_WHITELISTED_PROPERTY_LIST = ${CMAKE_WHITELISTED_PROPERTY_LIST}")
endfunction(print_whitelisted_properties)

function(print_target_properties tgt)
    if(NOT TARGET ${tgt})
        message("There is no target named '${tgt}'")
        return()
    endif()

    get_target_property(target_type ${tgt} TYPE)
    if(target_type STREQUAL "INTERFACE_LIBRARY")
        set(PROP_LIST ${CMAKE_WHITELISTED_PROPERTY_LIST})
    else()
        set(PROP_LIST ${CMAKE_PROPERTY_LIST})
    endif()

    foreach (prop ${PROP_LIST})
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
        # message ("Checking ${prop}")
        get_property(propval TARGET ${tgt} PROPERTY ${prop} SET)
        if (propval)
            get_target_property(propval ${tgt} ${prop})
            message ("${tgt} ${prop} = ${propval}")
        endif()
    endforeach(prop)
endfunction(print_target_properties)

print_target_properties(ga::ga)
print_target_properties(ga_src)
print_target_properties(ppidd)
print_target_properties(ppidd::ppidd)
