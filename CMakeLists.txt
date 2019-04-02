cmake_minimum_required(VERSION 3.9)

#
# libboink
#
project(boink VERSION 0.1 DESCRIPTION "streaming cDBG and dBG sketching algorithms")

# boink headers are listed in a cmake file.
include(${CMAKE_SOURCE_DIR}/cmake/boink_headers.cmake)
include(GNUInstallDirs)

# we want the conda includes if we're in a conda environment.
include_directories($ENV{CONDA_PREFIX}/include)

#
# Add our project's cmake dir the the module path. This gives us the
# Cppyy commands and targets.
#
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)
find_package(Cppyy)

#
# Make the default build us c++14 and "RELEASE" (-O3)
#
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE) 
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

#
# Add CQF as an object target: we'll include the resulting object directly
# in the boink shared library build
#
add_library(cqf OBJECT
            src/cqf/gqf.c 
)
target_include_directories(cqf PRIVATE include/cqf)

#
# Do the same object target method for smhasher.
#
add_library(smhasher OBJECT
            src/smhasher/MurmurHash3.cc
)
target_include_directories(smhasher PRIVATE include/smhasher)

#
# gfakluge has its own Makefile for produces its static library; run it.
#
add_custom_target(gfakluge
                  COMMAND make
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/third-party/gfakluge
)

#
# The boink shared library needs all the source files, the direct object
# targets, and the gfakluge static library.
#
add_library(boink SHARED
            src/boink/storage/partitioned_storage.cc
            src/boink/storage/storage.cc
            src/boink/storage/qfstorage.cc
            src/boink/storage/bytestorage.cc
            src/boink/storage/bitstorage.cc
            src/boink/storage/sparseppstorage.cc
            src/boink/storage/nibblestorage.cc
            src/boink/processors.cc
            src/boink/metrics.cc
            src/boink/events.cc
            src/boink/reporting/cdbg_writer_reporter.cc
            src/boink/reporting/ukhs_signature_reporter.cc
            src/boink/reporting/report_types.cc
            src/boink/reporting/cdbg_unitig_reporter.cc
            src/boink/reporting/cdbg_history_reporter.cc
            src/boink/reporting/reporters.cc
            src/boink/reporting/cdbg_component_reporter.cc
            src/boink/reporting/streaming_compactor_reporter.cc
            src/boink/hashing/hashshifter.cc
            src/boink/hashing/alphabets.cc
            src/boink/hashing/kmeriterator.cc
            src/boink/hashing/rollinghashshifter.cc
            src/boink/hashing/ukhs.cc
            src/boink/hashing/hashing_types.cc
            src/boink/ukhs_signature.cc
            src/boink/assembly.cc
            src/boink/dbg.cc
            src/boink/normalization/diginorm.cc
            src/boink/boink.cc
            src/boink/event_types.cc
            src/boink/cdbg/metrics.cc
            src/boink/cdbg/cdbg_types.cc
            src/boink/cdbg/cdbg.cc
            src/boink/cdbg/compactor.cc
            src/boink/cdbg/solid_compactor.cc
            src/boink/parsing/readers.cc
            src/boink/parsing/parsing.cc
            src/boink/minimizers.cc
            $<TARGET_OBJECTS:cqf>    
            $<TARGET_OBJECTS:smhasher>
)
set_target_properties(boink PROPERTIES
                      VERSION ${PROJECT_VERSION}
                      SOVERSION 1
)

#
# Prometheus has its own CMake files. We add it as a subdirectory, while also
# grabbing its include dirs (TODO: is there a way to pull them from the CMake
# exports before its installed?)
#
set(PROMETHEUS_DIR ${CMAKE_SOURCE_DIR}/third-party/prometheus-cpp)
add_subdirectory(${PROMETHEUS_DIR})
set(PROMETHEUS_INCLUDE_DIRS ${PROMETHEUS_DIR}/core/include ${PROMETHEUS_DIR}/pull/include ${PROMETHEUS_DIR}/push/include)

#
# Find zlib. This is a standard way of bringing it in.
#
find_package(ZLIB REQUIRED)
target_include_directories(boink PRIVATE ${ZLIB_INCLUDE_DIRS})

#
# Link the prometheus shared libs, zlib, and the generated gfakluge static
# lib to the boink shared lib.
#
target_link_libraries(boink
                      ${ZLIB_LIBRARIES} 
                      prometheus-cpp::pull
                      prometheus-cpp::push
                      prometheus-cpp::core
                      "${CMAKE_SOURCE_DIR}/third-party/gfakluge/libgfakluge.a"
)
# Make sure CMake knows that boink depends on the gfakluge custom command target
add_dependencies(boink gfakluge)

#
# Add all the include directories. Includes that are only used in boink's
# source files are set private, so that they are not installed.
#
target_include_directories(boink PUBLIC third-party/gfakluge/src/)
target_include_directories(boink PRIVATE third-party/)
target_include_directories(boink PRIVATE include/)
target_include_directories(boink PUBLIC include/boink)
target_include_directories(boink PUBLIC ${PROMETHEUS_INCLUDE_DIRS})

#
# Configure the boink pc.
#
configure_file(src/boink/boink.pc.in boink.pc @ONLY)

#
# libboink's install commands. Installs the libboink shared so,
# its headers, the prometheus headers, and any other public includes.
#
install(TARGETS boink
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/boink/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/boink
)

install(FILES ${CMAKE_BINARY_DIR}/boink.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig
)

install(FILES third-party/gfakluge/src/gfakluge.hpp third-party/gfakluge/src/pliib.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

#
# Set up the Cppyy bindings generation. This is a customized version defined
# in boink's cmake/ dir; it uses genreflex rather than calling rootcling directly.
# I did this because I couldn't get rootcling to properly include/exclude classes
# via the LinkDef header, and I wanted to be able to use the better syntax in
# the genreflex selection XML anyhow. Also, I think this is now the recommended /
# more modern way anyhow? Code was modified from the versions cppyy distributes.
#
cppyy_add_bindings(
    "boink" "${PROJECT_VERSION}" "Camille Scott" "cswel@ucdavis.edu"
   LICENSE "MIT"
   LANGUAGE_STANDARD "14"
   SELECTION_XML  ${CMAKE_CURRENT_SOURCE_DIR}/boink_select.xml
   INTERFACE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/include/boink/interface.hh
   HEADERS        ${BOINK_HEADERS}
   INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}
                  ${CMAKE_CURRENT_SOURCE_DIR}/include
                  ${CMAKE_CURRENT_SOURCE_DIR}/third-party
                  ${CMAKE_CURRENT_SOURCE_DIR}/third-party/gfakluge/src
                  ${ZLIB_INCLUDE_DIRS}
                  ${PROMETHEUS_INCLUDE_DIRS}
   LINK_LIBRARIES boink
)
