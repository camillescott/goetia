# boink/tests/test_dbg.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from boink.tests.utils import *
from boink.processors import FileConsumer

from boink.dbg import PdBG

@using_ksize(21)
def test_K(graph, ksize):
    assert graph.K == ksize


@using_ksize([21, 51, 101])
@presence_backends()
def test_presence(graph, ksize, random_sequence):
    # basic get/add test
    for kmer in kmers(random_sequence(), ksize):

        hashval = graph.hash(kmer)

        assert graph.query(kmer) == 0
        assert graph.query(hashval) == 0

        graph.insert(kmer)
        assert graph.query(kmer) == 1
        assert graph.query(hashval) == 1

        graph.insert(kmer)
        assert graph.query(kmer) == 1
        assert graph.query(hashval) == 1


@using_ksize([21, 51, 101])
@counting_backends()
def test_counting_presence(graph, ksize, random_sequence):
    # basic get/add test
    for kmer in kmers(random_sequence(), ksize):

        hashval = graph.hash(kmer)

        assert graph.query(kmer) == 0
        assert graph.query(hashval) == 0

        graph.insert(kmer)
        assert graph.query(kmer) == 1
        assert graph.query(hashval) == 1

        graph.insert(kmer)
        assert graph.query(kmer) == 2
        assert graph.query(hashval) == 2


@using_ksize([21, 51, 101])
@counting_backends()
def test_counting_count(graph, ksize, random_sequence):
    seq_kmers = list(kmers(random_sequence(), ksize))
    for iterations in range(10):
        for kmer in seq_kmers:
            hashval = graph.hash(kmer)
            assert graph.query(hashval) == iterations
            assert graph.query(kmer)    == iterations
            graph.insert(hashval)


@using_ksize([21, 51, 101])
@counting_backends()
def test_counting_count_add_sequence(graph, ksize, random_sequence):
    seq = random_sequence()
    seq_kmers = list(kmers(seq, ksize))
    for iterations in range(10):
        for kmer in seq_kmers:
            hashval = graph.hash(kmer)
            assert graph.query(hashval) == iterations
            assert graph.query(kmer)    == iterations
        graph.insert_sequence(seq)


@using_ksize([21,151])
@oxli_backends()
def test_n_occupied(graph, ksize):
    # basic get/add test
    kmer = 'G' * ksize

    assert graph.n_occupied == 0
    assert graph.n_unique == 0

    graph.insert(kmer)
    assert graph.n_occupied == 1
    assert graph.n_unique == 1

    graph.insert(kmer)
    # the CQF implementation we use can use more than one slot to represent
    # counts for a single kmer
    if not "QF" in graph.__class__.__name__:
        assert graph.n_occupied == 1
    else:
        assert graph.n_occupied == 2
    assert graph.n_unique == 1


@using_ksize([21,51,81])
def test_get_ksize(graph, ksize):
    assert graph.K == ksize


@using_ksize(5)
def test_hash(graph, ksize):
    # hashing of strings -> numbers.
    x = graph.hash("ATGGC")
    assert type(x) == int


def test_hash_bad_dna(graph, ksize):
    # hashing of bad dna -> succeeds w/o complaint
    with pytest.raises(ValueError):
        x = graph.hash("ATGYC")


def test_hash_bad_length(graph, ksize):
    # hashing of too long should ignore extra sequence
    test_kmer = 'A' * ksize
    assert graph.hash(test_kmer) == graph.hash(test_kmer + 'TTTT')


@using_ksize(5)
def test_add_hashval(graph, ksize):
    # test add(hashval)
    x = graph.hash("ATGGC")
    y = graph.insert(x)
    assert y

    z = graph.query(x)
    assert z == 1


@using_ksize(5)
def test_add_dna_kmer(graph, ksize):
    # test add(dna)
    x = graph.insert("ATGGC")
    assert x

    z = graph.query("ATGGC")
    assert z == 1


@using_ksize(5)
def test_get_hashval(graph, ksize):
    # test get(hashval)
    hashval = graph.hash("ATGGC")
    graph.insert(hashval)

    z = graph.query(hashval)
    assert z == 1


@using_ksize(5)
def test_get_dna_kmer(graph, ksize):
    # test get(dna)
    hashval = graph.hash("ATGGC")
    graph.insert(hashval)

    z = graph.query("ATGGC")
    assert z == 1


@using_ksize(5)
def test_get_bad_dna_kmer(graph, ksize):
    # test get(dna) with bad dna; should fail
    with pytest.raises(ValueError):
        graph.query("ATYGC")


@using_ksize(5)
def test_neighbors(graph, ksize):
    s = "ATGCCGATGCA"
    graph.insert_sequence(s)
    l, r = graph.neighbors(s[:5])
    assert len(l) == 1
    assert len(r) == 1


@using_ksize(5)
def test_add_sequence_and_report(graph, ksize):
    x = "ATGCCGATGCA"
    _, report = graph.insert_sequence_and_report(x)
    num_kmers = sum(report)
    assert num_kmers == len(x) - ksize + 1   # num k-mers consumed

    for start in range(len(x) - 6 + 1):
        assert graph.query(x[start:start + 6]) == 1


@using_ksize(5)
def test_add_sequence_bad_dna(graph):
    # while we don't specifically handle bad DNA, we should at least be
    # consistent...
    x = "ATGCCGNTGCA"
    with pytest.raises(ValueError):
        num_kmers = graph.insert_sequence(x)


@using_ksize(10)
def test_add_sequence_short(graph):
    # raise error on too short when consume is run
    x = "ATGCA"
    with pytest.raises(ValueError):
        graph.insert_sequence(x)


@using_ksize(6)
def test_get_kmer_counts(graph):
    graph.insert_sequence("AAAAAA")
    counts = graph.query_sequence("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] == 1

    graph.insert_sequence("AAAAAA")
    counts = graph.query_sequence("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] >= 1

    graph.insert_sequence("AAAAAT")
    counts = graph.query_sequence("AAAAAAT")
    print(counts)
    assert len(counts) == 2
    assert counts[0] >= 1
    assert counts[1] == 1


@using_ksize(5)
def test_consume_alias(graph_type, ksize):
    _, graph_type = graph_type
    g1 = graph_type(ksize, 1000, 4)
    g2 = graph_type(ksize, 1000, 4)
    assert g1.insert_sequence('A' * (ksize + 1)) == \
           g2.consume('A' * (ksize + 1))


@using_ksize(6)
def test_get_kmer_hashes(graph):
    hashes = list(graph.hashes("ACGTGCGT"))
    print(hashes)
    assert len(hashes) == 3
    assert hashes[0] == graph.hash("ACGTGC")
    assert hashes[1] == graph.hash("CGTGCG")
    assert hashes[2] == graph.hash("GTGCGT")


@using_ksize([21, 31, 41])
@using_length(1000)
def test_hashing_2(graph, linear_path, ksize):
    ''' Graph.hash uses a stand alone hasher for RollingHashShifters,
    Graph.hashes uses a KmerIterator; check that they give the same
    results.'''

    S = linear_path()

    for u, v in zip((graph.hash(kmer) for kmer in kmers(S, ksize)),
                    graph.hashes(S)):
        assert u == v


@using_ksize([21, 31, 41])
@using_length([50000, 500000])
@pytest.mark.benchmark(group='dbg-sequence')
@exact_backends()
def test_n_unique(graph, random_sequence, ksize, benchmark):
    sequence = random_sequence()
    kmer_set = set(kmers(sequence, ksize))
    benchmark(graph.insert_sequence, sequence)

    assert len(kmer_set) == graph.n_unique


@using_ksize([21, 31, 41])
@using_length([50000, 500000])
@pytest.mark.benchmark(group='dbg-sequence')
def test_get_counts(graph, random_sequence, ksize, benchmark):
    sequence = random_sequence()
    graph.insert_sequence(sequence)

    counts = benchmark(graph.query_sequence, sequence)
    assert all((count > 0 for count in counts))


@using_ksize([21, 31, 41])
@using_length([50000, 500000])
@pytest.mark.benchmark(group='dbg-sequence')
def test_pdbg_n_unique(random_sequence, ksize, benchmark):
    graph = PdBG(ksize, 7)
    sequence = random_sequence()
    kmer_set = set(kmers(sequence, ksize))
    benchmark(graph.insert_sequence, sequence)

    assert len(kmer_set) == graph.n_unique


@using_ksize([21, 31, 41])
@using_length([50000, 500000])
@pytest.mark.benchmark(group='dbg-sequence')
def test_pdbg_get_counts(random_sequence, ksize, benchmark):
    sequence = random_sequence()
    graph = PdBG(ksize, 7)
    graph.insert_sequence(sequence)

    counts = benchmark(graph.query_sequence, sequence)
    assert all((count > 0 for count in counts))
