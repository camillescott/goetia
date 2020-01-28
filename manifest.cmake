set(_headers
    include/boink/benchmarks/bench_storage.hh
    include/boink/boink.hh
    include/boink/cdbg/cdbg.hh
    include/boink/cdbg/cdbg_types.hh
    include/boink/cdbg/compactor.hh
    include/boink/cdbg/metrics.hh
    include/boink/cdbg/saturating_compactor.hh
    include/boink/cdbg/ucompactor.hh
    include/boink/cdbg/udbg.hh
    include/boink/cdbg/utagger.hh
    include/boink/dbg.hh
    include/boink/events.hh
    include/boink/event_types.hh
    include/boink/sequences/alphabets.hh
    include/boink/hashing/hash_combine.hh
    include/boink/hashing/canonical.hh
    include/boink/hashing/hashshifter.hh
    include/boink/hashing/hashextender.hh
    include/boink/hashing/kmeriterator.hh
    include/boink/hashing/kmer_span.hh
    include/boink/hashing/shifter_types.hh
    include/boink/hashing/rollinghash/characterhash.h
    include/boink/hashing/rollinghash/cyclichash.h
    include/boink/hashing/rollinghashshifter.hh
    include/boink/hashing/smhasher/MurmurHash3.h
    include/boink/hashing/unikmershifter.hh
    include/boink/hashing/ukhs.hh
    include/boink/interface.hh
    include/boink/is_detected.hh
    include/boink/meta.hh
    include/boink/metrics.hh
    include/boink/minimizers.hh
    include/boink/parsing/kseq.h
    include/boink/parsing/parsing.hh
    include/boink/parsing/readers.hh
    include/boink/pdbg.hh
    include/boink/processors.hh
    include/boink/reporting/reporters.hh
    include/boink/ring_span.hpp
    include/boink/signatures/sourmash/kmer_min_hash.hh
    include/boink/signatures/sourmash_signature.hh
    include/boink/signatures/ukhs_signature.hh
    include/boink/storage/bitstorage.hh
    include/boink/storage/bytestorage.hh
    include/boink/storage/cqf/gqf.h
    include/boink/storage/nibblestorage.hh
    include/boink/storage/partitioned_storage.hh
    include/boink/storage/qfstorage.hh
    include/boink/storage/sparsepp/serialize.hh
    include/boink/storage/sparseppstorage.hh
    include/boink/storage/storage.hh
    include/boink/storage/storage_types.hh
    include/boink/traversal.hh
    include/boink/utils/stringutils.h
)

set(_sources
    src/boink/pdbg.cc
    src/boink/storage/qfstorage.cc
    src/boink/storage/bytestorage.cc
    src/boink/storage/bitstorage.cc
    src/boink/storage/sparseppstorage.cc
    src/boink/storage/nibblestorage.cc
    src/boink/signatures/ukhs_signature.cc
    src/boink/signatures/sourmash_signature.cc
    src/boink/benchmarks/bench_storage.cc
    src/boink/events.cc
    src/boink/reporting/reporters.cc
    src/boink/hashing/hashshifter.cc
    src/boink/hashing/hashextender.cc
    src/boink/hashing/alphabets.cc
    src/boink/hashing/kmeriterator.cc
    src/boink/hashing/kmer_span.cc
    src/boink/hashing/rollinghashshifter.cc
    src/boink/hashing/unikmershifter.cc
    src/boink/hashing/smhasher/MurmurHash3.cc
    src/boink/hashing/ukhs.cc
    src/boink/hashing/canonical.cc
    src/boink/sequences/alphabets.cc
    src/boink/dbg.cc
    src/boink/traversal.cc
    src/boink/boink.cc
    src/boink/meta.cc
    src/boink/cdbg/metrics.cc
    src/boink/cdbg/cdbg.cc
    src/boink/cdbg/compactor.cc
    src/boink/cdbg/ucompactor.cc
    src/boink/cdbg/utagger.cc
    src/boink/cdbg/udbg.cc
    src/boink/cdbg/saturating_compactor.cc
    src/boink/parsing/readers.cc
    src/boink/parsing/parsing.cc
    src/boink/minimizers.cc
    src/boink/storage/cqf/gqf.c
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
    set(BOINK_INCLUDE_ROOT $ENV{CONDA_PREFIX})
else()
    set(BOINK_INCLUDE_ROOT ${CMAKE_SOURCE_DIR})
endif()

foreach (path ${_headers})
    list(APPEND LIB_HEADERS ${BOINK_INCLUDE_ROOT}/${path})
endforeach(path)

foreach (path ${_sources})
    list(APPEND LIB_SOURCES ${CMAKE_SOURCE_DIR}/${path})
endforeach(path)

foreach (path ${_data})
    list(APPEND LIB_DATA ${CMAKE_SOURCE_DIR}/py/data/${path})
endforeach(path)
