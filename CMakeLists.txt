cmake_minimum_required(VERSION 3.9)
project(boink VERSION 0.1 DESCRIPTION "streaming cDBG and dBG sketching algorithms")
include(GNUInstallDirs)
include_directories($ENV{CONDA_PREFIX}/include)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


add_library(cqf OBJECT
            src/cqf/gqf.c 
)
target_include_directories(cqf PRIVATE include/cqf)

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
)

set_target_properties(boink PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
)

add_subdirectory(third-party/prometheus-cpp)

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


target_include_directories(boink PUBLIC third-party/sparsepp/)
target_include_directories(boink PUBLIC third-party/rollinghash/)
target_include_directories(boink PUBLIC third-party/gfakluge/src/)

target_include_directories(boink PRIVATE third-party/)
target_include_directories(boink PRIVATE third-party/gfakluge/)

configure_file(src/boink/boink.pc.in boink.pc @ONLY)
target_include_directories(boink PRIVATE include/)
target_include_directories(boink PUBLIC include/boink)

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
