cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0063 NEW)

#
# libgoetia
#

file(STRINGS "${CMAKE_SOURCE_DIR}/include/goetia/VERSION" sharedlib_version)

message(STATUS "version: ${sharedlib_version}")

include(CheckCXXCompilerFlag)
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
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

project(goetia
        LANGUAGES CXX  C
        VERSION ${sharedlib_version}
        DESCRIPTION "streaming cDBG and dBG sketching algorithms"
        HOMEPAGE_URL "https://github.com/camillescott/goetia"
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

# Find sourmash
#
find_library(LIBSOURMASH sourmash)
message(STATUS ${LIBSOURMASH})

#
# The goetia shared library needs all the source files, the direct object
# targets, and the gfakluge static library.
#
add_library(goetia SHARED
            ${LIB_SOURCES}
)
set_target_properties(goetia PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(goetia PROPERTIES CXX_VISIBILITY_PRESET default)
set_property(TARGET goetia PROPERTY LINK_WHAT_YOU_USE TRUE)
set_property(TARGET goetia PROPERTY VISIBILITY_INLINES_HIDDEN 0)
set_property(TARGET goetia PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
set_target_properties(goetia PROPERTIES
                      VERSION ${PROJECT_VERSION}
                      SOVERSION 1
)

#set_target_properties(goetia PROPERTIES PUBLIC_HEADER ${LIB_HEADERS})
target_include_directories(goetia
                           PUBLIC
                           ${CMAKE_SOURCE_DIR}/include
)
target_include_directories(goetia PRIVATE ${ZLIB_INCLUDE_DIRS})
target_link_libraries(goetia
                      Threads::Threads
                      ${ZLIB_LIBRARIES} 
                      ${LIBSOURMASH}
)


#
# Add all the include directories. Includes that are only used in goetia's
# source files are set private, so that they are not installed.
#
target_include_directories(goetia PRIVATE third-party/)
target_include_directories(goetia PRIVATE include/)
target_include_directories(goetia PRIVATE include/goetia/sketches/sketch/)
target_include_directories(goetia PRIVATE include/goetia/sketches/sketch/include)

#
# Configure the goetia pc.
#
configure_file(src/goetia/goetia.pc.in goetia.pc @ONLY)


#
# Benchmark exes
#
add_executable(do_bench_storage ${CMAKE_SOURCE_DIR}/src/goetia/benchmarks/do_bench_storage.cc)
target_link_libraries(do_bench_storage goetia)

add_executable(test_hashing EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/src/goetia/benchmarks/benchmark_hashing.cc)
target_link_libraries(test_hashing goetia)

#
# Set up the Cppyy bindings generation. This is a customized version defined
# in goetia's cmake/ dir; it uses genreflex rather than calling rootcling directly.
# I did this because I couldn't get rootcling to properly include/exclude classes
# via the LinkDef header, and I wanted to be able to use the better syntax in
# the genreflex selection XML anyhow. Also, I think this is now the recommended /
# more modern way anyhow? Code was modified from the versions cppyy distributes.
#
cppyy_add_bindings(
   "goetia" 
   LANGUAGE_STANDARD "17"
   SELECTION_XML  ${GOETIA_INCLUDE_ROOT}/include/goetia/interface.xml
   INTERFACE_FILE ${GOETIA_INCLUDE_ROOT}/include/goetia/interface.hh
   HEADERS        ${LIB_HEADERS}
   INTERFACE_HEADERS ${API_HEADERS}
   INCLUDE_DIRS   ${GOETIA_INCLUDE_ROOT}/include
                  ${GOETIA_INCLUDE_ROOT}/include/goetia/sketches/sketch/
                  ${GOETIA_INCLUDE_ROOT}/include/goetia/sketches/sketch/include/
                  ${CMAKE_CURRENT_SOURCE_DIR}/third-party
                  ${ZLIB_INCLUDE_DIRS}
   LINK_LIBRARIES goetia
)

#
# libgoetia's install commands. Installs the libgoetia shared so,
# its headers, and any other public includes.
#
install(TARGETS goetia
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(TARGETS ${CPPYY_LIB_TARGET}
        LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/goetia)

install(FILES  ${CPPYY_ROOTMAP} ${CPPYY_EXTRA_MAP} ${CPPYY_PCM}
        DESTINATION ${CMAKE_SOURCE_DIR}/goetia)

install(TARGETS do_bench_storage
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR/benchmarks})
#install(DIRECTORY include/goetia/
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/goetia
#)

install(FILES ${CMAKE_BINARY_DIR}/goetia.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig
)

