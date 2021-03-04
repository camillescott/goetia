set(_headers
    include/goetia/benchmarks/bench_storage.hh
    include/goetia/goetia.hh
    include/goetia/cdbg/cdbg.hh
    include/goetia/cdbg/cdbg_types.hh
    include/goetia/cdbg/compactor.hh
    include/goetia/cdbg/metrics.hh
    include/goetia/cdbg/saturating_compactor.hh
    include/goetia/cdbg/ucompactor.hh
    include/goetia/cdbg/udbg.hh
    include/goetia/cdbg/utagger.hh
    include/goetia/dbg.hh
    include/goetia/sequences/alphabets.hh
    include/goetia/hashing/hash_combine.hh
    include/goetia/hashing/canonical.hh
    include/goetia/hashing/hashshifter.hh
    include/goetia/hashing/hashextender.hh
    include/goetia/hashing/kmeriterator.hh
    include/goetia/hashing/kmer_span.hh
    include/goetia/hashing/shifter_types.hh
    include/goetia/hashing/rollinghash/characterhash.h
    include/goetia/hashing/rollinghash/cyclichash.h
    include/goetia/hashing/rollinghashshifter.hh
    include/goetia/hashing/smhasher/MurmurHash3.h
    include/goetia/hashing/unikmershifter.hh
    include/goetia/hashing/ukhs.hh
    include/goetia/interface.hh
    include/goetia/is_detected.hh
    include/goetia/meta.hh
    include/goetia/metrics.hh
    include/goetia/minimizers.hh
    include/goetia/parsing/kseq.h
    include/goetia/parsing/parsing.hh
    include/goetia/parsing/readers.hh
    include/goetia/pdbg.hh
    include/goetia/processors.hh
    include/goetia/ring_span.hpp
    include/goetia/solidifier.hh
    include/goetia/signatures/sourmash/sourmash.hpp
    include/goetia/signatures/sourmash/sourmash.h
    include/goetia/signatures/sourmash_signature.hh
    include/goetia/signatures/ukhs_signature.hh
    include/goetia/storage/bitstorage.hh
    include/goetia/storage/bytestorage.hh
    include/goetia/storage/cqf/gqf.h
    include/goetia/storage/nibblestorage.hh
    include/goetia/storage/partitioned_storage.hh
    include/goetia/storage/qfstorage.hh
    include/goetia/storage/sparsepp/serialize.hh
    include/goetia/storage/sparseppstorage.hh
    include/goetia/storage/storage.hh
    include/goetia/storage/storage_types.hh
    include/goetia/traversal.hh
    include/goetia/utils/stringutils.h
)

set(_sources
    src/goetia/pdbg.cc
    src/goetia/storage/qfstorage.cc
    src/goetia/storage/bytestorage.cc
    src/goetia/storage/bitstorage.cc
    src/goetia/storage/sparseppstorage.cc
    src/goetia/storage/nibblestorage.cc
    src/goetia/signatures/ukhs_signature.cc
    src/goetia/signatures/sourmash_signature.cc
    src/goetia/benchmarks/bench_storage.cc
    src/goetia/hashing/hashshifter.cc
    src/goetia/hashing/hashextender.cc
    src/goetia/hashing/kmeriterator.cc
    src/goetia/hashing/kmer_span.cc
    src/goetia/hashing/rollinghashshifter.cc
    src/goetia/hashing/unikmershifter.cc
    src/goetia/hashing/smhasher/MurmurHash3.cc
    src/goetia/hashing/ukhs.cc
    src/goetia/hashing/canonical.cc
    src/goetia/sequences/alphabets.cc
    src/goetia/dbg.cc
    src/goetia/traversal.cc
    src/goetia/solidifier.cc
    src/goetia/goetia.cc
    src/goetia/meta.cc
    src/goetia/processors.cc
    src/goetia/cdbg/metrics.cc
    src/goetia/cdbg/cdbg.cc
    src/goetia/cdbg/compactor.cc
    src/goetia/cdbg/ucompactor.cc
    src/goetia/cdbg/utagger.cc
    src/goetia/cdbg/udbg.cc
    src/goetia/cdbg/saturating_compactor.cc
    src/goetia/parsing/readers.cc
    src/goetia/parsing/parsing.cc
    src/goetia/minimizers.cc
    src/goetia/storage/cqf/gqf.c
)

set(_data
    res_10_100_4_0.txt.gz res_10_110_4_0.txt.gz res_10_120_4_0.txt.gz
    res_10_130_4_0.txt.gz res_10_140_4_0.txt.gz res_10_150_4_0.txt.gz
    res_10_160_4_0.txt.gz res_10_170_4_0.txt.gz res_10_180_4_0.txt.gz
    res_10_190_4_0.txt.gz res_10_200_4_0.txt.gz res_10_20_4_0.txt.gz
    res_10_30_4_0.txt.gz res_10_40_4_0.txt.gz res_10_50_4_0.txt.gz
    res_10_60_4_0.txt.gz res_10_70_4_0.txt.gz res_10_80_4_0.txt.gz
    res_10_90_4_0.txt.gz res_7_100_4_0.txt.gz res_7_110_4_0.txt.gz
    res_7_120_4_0.txt.gz res_7_130_4_0.txt.gz res_7_140_4_0.txt.gz
    res_7_150_4_0.txt.gz res_7_160_4_0.txt.gz res_7_170_4_0.txt.gz
    res_7_180_4_0.txt.gz res_7_190_4_0.txt.gz res_7_200_4_0.txt.gz
    res_7_20_4_0.txt.gz res_7_30_4_0.txt.gz res_7_40_4_0.txt.gz
    res_7_50_4_0.txt.gz res_7_60_4_0.txt.gz res_7_70_4_0.txt.gz
    res_7_80_4_0.txt.gz res_7_90_4_0.txt.gz res_8_100_4_0.txt.gz
    res_8_110_4_0.txt.gz res_8_120_4_0.txt.gz res_8_130_4_0.txt.gz
    res_8_140_4_0.txt.gz res_8_150_4_0.txt.gz res_8_160_4_0.txt.gz
    res_8_170_4_0.txt.gz res_8_180_4_0.txt.gz res_8_190_4_0.txt.gz
    res_8_200_4_0.txt.gz res_8_20_4_0.txt.gz res_8_30_4_0.txt.gz
    res_8_40_4_0.txt.gz res_8_50_4_0.txt.gz res_8_60_4_0.txt.gz
    res_8_70_4_0.txt.gz res_8_80_4_0.txt.gz res_8_90_4_0.txt.gz
    res_9_100_4_0.txt.gz res_9_110_4_0.txt.gz res_9_120_4_0.txt.gz
    res_9_130_4_0.txt.gz res_9_140_4_0.txt.gz res_9_150_4_0.txt.gz
    res_9_160_4_0.txt.gz res_9_170_4_0.txt.gz res_9_180_4_0.txt.gz
    res_9_190_4_0.txt.gz res_9_200_4_0.txt.gz res_9_20_4_0.txt.gz
    res_9_30_4_0.txt.gz res_9_40_4_0.txt.gz res_9_50_4_0.txt.gz
    res_9_60_4_0.txt.gz res_9_70_4_0.txt.gz res_9_80_4_0.txt.gz
    res_9_90_4_0.txt.gz
)

if(DEFINED ENV{CONDA_BUILD_DEPLOY})
    message(STATUS "Building a conda deployment, use installed headers.")
    set(GOETIA_INCLUDE_ROOT $ENV{CONDA_PREFIX})
else()
    set(GOETIA_INCLUDE_ROOT ${CMAKE_SOURCE_DIR})
endif()

foreach (path ${_headers})
    list(APPEND LIB_HEADERS ${GOETIA_INCLUDE_ROOT}/${path})
endforeach(path)

foreach (path ${_sources})
    list(APPEND LIB_SOURCES ${CMAKE_SOURCE_DIR}/${path})
endforeach(path)

foreach (path ${_data})
    list(APPEND LIB_DATA ${CMAKE_SOURCE_DIR}/py/data/${path})
endforeach(path)
