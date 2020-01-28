cmake_minimum_required(VERSION 3.9)
cmake_policy(SET CMP0063 NEW)

#
# libboink
#

if(DEFINED ENV{BOINK_VERSION})
    set(sharedlib_version "$ENV{BOINK_VERSION}")
    set(python_pkg_version "$ENV{BOINK_VERSION}")
else()
    execute_process(COMMAND python ../version.py --cmake
                    OUTPUT_VARIABLE sharedlib_version
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    )
    execute_process(COMMAND python ../version.py
                    OUTPUT_VARIABLE python_pkg_version
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    )
endif()
message(STATUS "version: ${python_pkg_version}")

project(boink
        LANGUAGES CXX  C
        VERSION ${sharedlib_version}
        DESCRIPTION "streaming cDBG and dBG sketching algorithms"
        HOMEPAGE_URL "https://github.com/camillescott/boink"
        )

if(DEFINED ENV{CONDA_PREFIX})
    message(STATUS "Building in a conda environment.")
    set(CONDA_ACTIVE TRUE)
    set(CMAKE_INSTALL_PREFIX "$ENV{CONDA_PREFIX}")
    set(CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
    set(CMAKE_INCLUDE_PATH "$ENV{CONDA_PREFIX}/include")
    set(CMAKE_LIBRARY_PATH "$ENV{CONDA_PREFIX}/lib")
endif()

include(${CMAKE_SOURCE_DIR}/manifest.cmake)
include(GNUInstallDirs)

set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads REQUIRED)


#
# Add our project's cmake dir the the module path. This gives us the
# Cppyy commands and targets.
#
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)
find_package(Cppyy)

# Find zlib. This is a standard way of bringing it in.
#
find_package(ZLIB REQUIRED)

#
# Make the default build use c++17 and "RELEASE" (-O3)
#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 0)
set(CMAKE_LINK_WHAT_YOU_USE TRUE)
if(NOT CMAKE_BUILD_TYPE) 
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release"
        FORCE
    )
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")

#
# The boink shared library needs all the source files, the direct object
# targets, and the gfakluge static library.
#
add_library(boink SHARED
            ${LIB_SOURCES}
)
set_target_properties(boink PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(boink PROPERTIES CXX_VISIBILITY_PRESET default)
set_property(TARGET boink PROPERTY LINK_WHAT_YOU_USE TRUE)
set_property(TARGET boink PROPERTY VISIBILITY_INLINES_HIDDEN 0)
set_target_properties(boink PROPERTIES
                      VERSION ${PROJECT_VERSION}
                      SOVERSION 1
)

#set_target_properties(boink PROPERTIES PUBLIC_HEADER ${LIB_HEADERS})
target_include_directories(boink
                           PUBLIC
                           ${CMAKE_SOURCE_DIR}/include
)
target_include_directories(boink PRIVATE ${ZLIB_INCLUDE_DIRS})
target_link_libraries(boink
                      Threads::Threads
                      ${ZLIB_LIBRARIES} 
)


#
# Add all the include directories. Includes that are only used in boink's
# source files are set private, so that they are not installed.
#
target_include_directories(boink PRIVATE third-party/)
target_include_directories(boink PRIVATE include/)

#
# Configure the boink pc.
#
configure_file(src/boink/boink.pc.in boink.pc @ONLY)


#
# Benchmark exes
#
add_executable(do_bench_storage EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/src/boink/benchmarks/do_bench_storage.cc)
target_link_libraries(do_bench_storage boink)

add_executable(test_hashing EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/src/boink/benchmarks/benchmark_hashing.cc)
target_link_libraries(test_hashing boink)

#
# Set up the Cppyy bindings generation. This is a customized version defined
# in boink's cmake/ dir; it uses genreflex rather than calling rootcling directly.
# I did this because I couldn't get rootcling to properly include/exclude classes
# via the LinkDef header, and I wanted to be able to use the better syntax in
# the genreflex selection XML anyhow. Also, I think this is now the recommended /
# more modern way anyhow? Code was modified from the versions cppyy distributes.
#
cppyy_add_bindings(
    "boink" ${python_pkg_version} "Camille Scott" "cswel@ucdavis.edu"
   LICENSE "MIT"
   LANGUAGE_STANDARD "17"
   SELECTION_XML  ${CMAKE_SOURCE_DIR}/boink_select.xml
   INTERFACE_FILE ${BOINK_INCLUDE_ROOT}/include/boink/interface.hh
   PKG_SRC_DIR    ${CMAKE_SOURCE_DIR}/boink
   TESTS_DIR      ${CMAKE_SOURCE_DIR}/tests
   HEADERS        ${LIB_HEADERS}
   INCLUDE_DIRS   ${BOINK_INCLUDE_ROOT}/include
                  ${CMAKE_CURRENT_SOURCE_DIR}/third-party
                  ${ZLIB_INCLUDE_DIRS}
   LINK_LIBRARIES boink
)

#
# libboink's install commands. Installs the libboink shared so,
# its headers, the prometheus headers, and any other public includes.
#
install(TARGETS boink
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

#install(DIRECTORY include/boink/
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/boink
#)

install(FILES ${CMAKE_BINARY_DIR}/boink.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig
)

install(CODE "execute_process(COMMAND python -m pip install --no-deps -U ${PY_WHEEL_FILE})")
