set(_headers
    include/boink/assembly.hh
    include/boink/boink.hh
    include/boink/cdbg/cdbg.hh
    include/boink/cdbg/cdbg_types.hh
    include/boink/cdbg/compactor.hh
    include/boink/cdbg/metrics.hh
    include/boink/dbg.hh
    include/boink/events.hh
    include/boink/event_types.hh
    include/boink/hashing/alphabets.hh
    include/boink/hashing/exceptions.hh
    include/boink/hashing/hashing_types.hh
    include/boink/hashing/hashshifter.hh
    include/boink/hashing/kmeriterator.hh
    include/boink/hashing/rollinghashshifter.hh
    include/boink/hashing/ukhshashshifter.hh
    include/boink/hashing/smhasher/MurmurHash3.h
    include/boink/hashing/ukhs.hh
    include/boink/interface.hh
    include/boink/kmers/kmerclient.hh
    include/boink/metrics.hh
    include/boink/minimizers.hh
    include/boink/parsing/parsing.hh
    include/boink/parsing/readers.hh
    include/boink/pdbg.hh
    include/boink/processors.hh
    include/boink/reporting/reporters.hh
    include/boink/storage/bitstorage.hh
    include/boink/storage/bytestorage.hh
    include/boink/storage/nibblestorage.hh
    include/boink/storage/partitioned_storage.hh
    include/boink/storage/qfstorage.hh
    include/boink/storage/sparseppstorage.hh
    include/boink/storage/storage.hh
    include/boink/storage/storage_types.hh
    include/boink/storage/cqf/gqf.h
    include/boink/signatures/ukhs_signature.hh
    include/boink/signatures/sourmash/kmer_min_hash.hh
)

set(_sources
    src/boink/storage/partitioned_storage.cc
    src/boink/storage/storage.cc
    src/boink/storage/qfstorage.cc
    src/boink/storage/bytestorage.cc
    src/boink/storage/bitstorage.cc
    src/boink/storage/sparseppstorage.cc
    src/boink/storage/nibblestorage.cc
    src/boink/storage/cqf/gqf.c
    src/boink/processors.cc
    src/boink/metrics.cc
    src/boink/events.cc
    src/boink/reporting/reporters.cc
    src/boink/hashing/smhasher/MurmurHash3.cc
    src/boink/hashing/hashshifter.cc
    src/boink/hashing/alphabets.cc
    src/boink/hashing/kmeriterator.cc
    src/boink/hashing/rollinghashshifter.cc
    src/boink/hashing/ukhs.cc
    src/boink/hashing/hashing_types.cc
    src/boink/signatures/ukhs_signature.cc
    src/boink/assembly.cc
    src/boink/dbg.cc
    src/boink/pdbg.cc
    src/boink/boink.cc
    src/boink/event_types.cc
    src/boink/cdbg/metrics.cc
    src/boink/cdbg/cdbg_types.cc
    src/boink/cdbg/cdbg.cc
    src/boink/cdbg/compactor.cc
    src/boink/parsing/readers.cc
    src/boink/parsing/parsing.cc
    src/boink/minimizers.cc
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

foreach (path ${_headers})
    list(APPEND LIB_HEADERS ${CMAKE_SOURCE_DIR}/${path})
endforeach(path)

foreach (path ${_sources})
    list(APPEND LIB_SOURCES ${CMAKE_SOURCE_DIR}/${path})
endforeach(path)

foreach (path ${_data})
    list(APPEND LIB_DATA ${CMAKE_SOURCE_DIR}/py/data/${path})
endforeach(path)
