cmake_minimum_required(VERSION 3.9)
project(boink VERSION 0.1 DESCRIPTION "streaming cDBG and dBG sketching algorithms")
include(boink_headers.cmake)
include(GNUInstallDirs)
include_directories($ENV{CONDA_PREFIX}/include)

execute_process(COMMAND cling-config --cmake OUTPUT_VARIABLE CPYY_MODULE_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
message("CPYY_MODULE_PATH: " ${CPYY_MODULE_PATH})
#list(INSERT CMAKE_MODULE_PATH 0 ${CPYY_MODULE_PATH})
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)
find_package(Cppyy)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


add_library(cqf OBJECT
            src/cqf/gqf.c 
)
target_include_directories(cqf PRIVATE include/cqf)


add_library(smhasher OBJECT
            src/smhasher/MurmurHash3.cc
)
target_include_directories(smhasher PRIVATE include/smhasher)


add_custom_target(gfakluge
                  COMMAND make
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/third-party/gfakluge
)


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

set(PROMETHEUS_DIR ${CMAKE_SOURCE_DIR}/third-party/prometheus-cpp)
add_subdirectory(${PROMETHEUS_DIR})
set(PROMETHEUS_INCLUDE_DIRS ${PROMETHEUS_DIR}/core/include ${PROMETHEUS_DIR}/pull/include ${PROMETHEUS_DIR}/push/include)

find_package(ZLIB REQUIRED)
target_include_directories(boink PRIVATE ${ZLIB_INCLUDE_DIRS})

target_link_libraries(boink
                      ${ZLIB_LIBRARIES} 
                      prometheus-cpp::pull
                      prometheus-cpp::push
                      prometheus-cpp::core
                      "${CMAKE_SOURCE_DIR}/third-party/gfakluge/libgfakluge.a"
)
add_dependencies(boink gfakluge)


target_include_directories(boink PUBLIC third-party/gfakluge/src/)
target_include_directories(boink PRIVATE third-party/)
target_include_directories(boink PRIVATE include/)
target_include_directories(boink PUBLIC include/boink)
target_include_directories(boink PUBLIC ${PROMETHEUS_INCLUDE_DIRS})


configure_file(src/boink/boink.pc.in boink.pc @ONLY)

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

#cppyy_add_bindings(
#   "boink" "0.1" "Camille Scott" "cswel@ucdavis.edu"
#   LICENSE "MIT"
#   LANGUAGE_STANDARD "14"
#   GENERATE_OPTIONS "-D__PIC__;-Wno-macro-redefined"
#   LINKDEFS boink_linkdef.h
#   INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
#                ${CMAKE_CURRENT_SOURCE_DIR}/include
#                ${CMAKE_CURRENT_SOURCE_DIR}/third-party
#                ${CMAKE_CURRENT_SOURCE_DIR}/third-party/gfakluge/src
#                ${ZLIB_INCLUDE_DIRS}
#                ${PROMETHEUS_INCLUDE_DIRS}
#
#   LINK_LIBRARIES boink
#   H_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include/boink
#   H_FILES ${BOINK_HEADER_NAMES}
#)


cppyy_add_bindings(
   "boink" "0.1" "Camille Scott" "cswel@ucdavis.edu"
   LICENSE "MIT"
   LANGUAGE_STANDARD "14"
   SELECTION_XML boink_select.xml
   INTERFACE_FILE include/boink/interface.hh
   INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/include
                ${CMAKE_CURRENT_SOURCE_DIR}/third-party
                ${CMAKE_CURRENT_SOURCE_DIR}/third-party/gfakluge/src
                ${ZLIB_INCLUDE_DIRS}
                ${PROMETHEUS_INCLUDE_DIRS}
   LINK_LIBRARIES boink
)
