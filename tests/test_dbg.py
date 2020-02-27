# goetia/tests/test_dbg.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from cppyy.gbl import std

from .utils import *
from goetia.hashing import FwdLemireShifter, UKHS
from goetia.storage import count_t
import pytest


@using(ksize=[21, 51, 101])
@presence_backends()
def test_presence_insert(graph, ksize, random_sequence):
    # basic get/add test
    for kmer in kmers(random_sequence(), ksize):

        hashval = graph.hash(kmer)

        assert graph.query(kmer) == 0
        assert graph.query(hashval) == 0

        assert graph.insert(kmer) == True
        assert graph.query(kmer) == 1
        assert graph.query(hashval) == 1

        assert graph.insert(kmer) == False
        assert graph.query(kmer) == 1
        assert graph.query(hashval) == 1


@using(ksize=[21, 51, 101])
@counting_backends()
def test_counting_insert(graph, ksize, random_sequence):
    # basic get/add test
    for kmer in kmers(random_sequence(), ksize):

        hashval = graph.hash(kmer)

        assert graph.query(kmer) == 0
        assert graph.query(hashval) == 0

        assert graph.insert(kmer) == True
        assert graph.query(kmer) == 1
        assert graph.query(hashval) == 1

        assert graph.insert(kmer) == False
        assert graph.query(kmer) == 2
        assert graph.query(hashval) == 2


@using(ksize=[21, 51, 101])
@counting_backends()
def test_counting_insert_and_query(graph, ksize, random_sequence):
    # basic get/add test
    for kmer in kmers(random_sequence(), ksize):

        hashval = graph.hash(kmer)
        assert graph.insert_and_query(kmer) == 1
        assert graph.insert_and_query(hashval) == 2


@using(ksize=[21, 51, 101])
@presence_backends()
def test_presence_insert_and_query(graph, ksize, random_sequence):
    # basic get/add test
    for kmer in kmers(random_sequence(), ksize):

        hashval = graph.hash(kmer)
        assert graph.insert_and_query(kmer) == 1
        assert graph.insert_and_query(hashval) == 1


@using(ksize=[21, 51, 101])
@counting_backends()
def test_counting_query(graph, ksize, random_sequence):
    seq_kmers = list(kmers(random_sequence(), ksize))
    for iterations in range(10):
        for kmer in seq_kmers:
            hashval = graph.hash(kmer)
            assert graph.query(hashval) == iterations
            assert graph.query(kmer)    == iterations
            graph.insert(hashval)


@using(ksize=[21, 51, 101])
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


@using(ksize=[21, 51, 101])
def test_get_ksize(graph, ksize):
    assert graph.K == ksize


@using(ksize=21)
def test_hash_type(graph, ksize, hasher_type):
    hasher_t, _ = hasher_type
    # hashing of strings -> numbers.
    x = graph.hash("A" * 21)
    assert issubclass(type(graph).hash_type, type(x))


@pytest.mark.xfail
def test_hash_bad_dna(graph, ksize):
    # hashing of bad dna -> succeeds w/o complaint
    # TODO: figure out cppyy exception conversion
    with pytest.raises(Exception):
        x = graph.hash("Y" * ksize)


def test_hash_bad_length(graph, ksize):
    # hashing of too long should ignore extra sequence
    test_kmer = 'A' * ksize
    assert graph.hash(test_kmer) == graph.hash(test_kmer + 'TTTT')


@using(ksize=21)
def test_add_hashval(graph, ksize):
    # test add(hashval)
    x = graph.hash("ATC" * 7)
    y = graph.insert(x)
    assert y

    z = graph.query(x)
    assert z == 1


@using(ksize=21)
def test_add_dna_kmer(graph, ksize):
    # test add(dna)
    x = graph.insert("ATC" * 7)
    assert x

    z = graph.query("ATC" * 7)
    assert z == 1


@using(ksize=21)
def test_get_hashval(graph, ksize):
    # test get(hashval)
    hashval = graph.hash("ATC" * 7)
    graph.insert(hashval)

    z = graph.query(hashval)
    assert z == 1


@using(ksize=21)
def test_get_dna_kmer(graph):
    # test get(dna)
    kmer = "ATC" * 7
    hashval = graph.hash(kmer)
    graph.insert(hashval)

    z = graph.query(kmer)
    assert z == 1


@using(ksize=21)
@pytest.mark.xfail
def test_get_bad_dna_kmer(graph, ksize):
    # test get(dna) with bad dna; should fail
    #TODO: figure out cppyy exception foo
    with pytest.raises(TypeError):
        graph.query("Y" * ksize)


@using(ksize=21, length=23)
def test_out_neighbors(graph, ksize, linear_path):
    s = linear_path()
    graph.insert_sequence(s)
    n = graph.out_neighbors(s[1:22])
    assert len(n) == 1


@using(ksize=21, length=23)
def test_neighbors(graph, ksize, linear_path):
    s = linear_path()
    graph.insert_sequence(s)
    l, r = graph.neighbors(s[1:22])
    assert len(l) == 1
    assert len(r) == 1


@using(ksize=21, length=50)
def test_insert_sequence_overload(graph, ksize, length, linear_path):
    x = linear_path()
    hashes = std.vector[type(graph).hash_type]()
    report = std.vector[count_t]()
    n_consumed = graph.insert_sequence(x, hashes, report)
    num_kmers = sum(report)
    assert num_kmers == len(x) - ksize + 1   # num k-mers consumed

    for start in range(len(x) - ksize + 1):
        assert graph.query(x[start:start + ksize]) == 1


@using(ksize=21, length=30)
@pytest.mark.xfail
def test_insert_sequence_bad_dna(graph, linear_path):
    # while we don't specifically handle bad DNA, we should at least be
    # consistent...
    x = linear_path() + 'X'
    with pytest.raises(Exception):
        num_kmers = graph.insert_sequence(x)


@using(ksize=21)
def test_add_sequence_short(graph):
    # raise error on too short when consume is run
    x = "ATGCA"
    with pytest.raises(Exception):
        graph.insert_sequence(x)


@using(ksize=21)
def test_get_kmer_counts(graph, ksize):
    graph.insert_sequence("A" * ksize)
    counts = graph.query_sequence("A" * ksize)
    print(counts)
    assert len(counts) == 1
    assert counts[0] == 1

    graph.insert_sequence("A" * ksize)
    counts = graph.query_sequence("A" * ksize)
    print(counts)
    assert len(counts) == 1
    assert counts[0] >= 1

    graph.insert_sequence(("A" * ksize) + 'T')
    counts = graph.query_sequence(("A" * ksize) + 'T')
    print(counts)
    assert len(counts) == 2
    assert counts[0] >= 1
    assert counts[1] == 1


@using(ksize=21, length=23)
def test_get_kmer_hashes(graph, ksize, length, linear_path):
    s = linear_path()
    hashes = list(graph.hashes(s))
    print(hashes)
    assert len(hashes) == 3
    assert hashes[0] == graph.hash(s[0:ksize])
    assert hashes[1] == graph.hash(s[1:1+ksize])
    assert hashes[2] == graph.hash(s[2:2+ksize])


@using(ksize=[21, 31, 41], length=1000)
def test_hashing_2(graph, linear_path, ksize):
    ''' Graph.hash uses a stand alone hasher for LemireShifterPolicys,
    Graph.hashes uses a KmerIterator; check that they give the same
    results.'''

    S = linear_path()

    for u, v in zip((graph.hash(kmer) for kmer in kmers(S, ksize)),
                    graph.hashes(S)):
        assert u == v


@using(ksize=[21, 31, 41], length=[50000, 500000])
@pytest.mark.benchmark(group='dbg-sequence')
@exact_backends()
def test_n_unique(graph, random_sequence, ksize, benchmark):
    sequence = random_sequence()
    kmer_set = set((graph.hash(kmer) for kmer in kmers(sequence, ksize)))
    benchmark(graph.insert_sequence, sequence)

    assert len(kmer_set) == graph.n_unique()


@using(ksize=[21, 31, 41], length=[50000, 500000])
@pytest.mark.benchmark(group='dbg-sequence')
def test_get_counts(graph, random_sequence, ksize, benchmark):
    sequence = random_sequence()
    graph.insert_sequence(sequence)

    counts = benchmark(graph.query_sequence, sequence)
    assert all((count > 0 for count in counts))


@using(ksize=[21, 51, 101])
def test_clone(graph, ksize, random_sequence):
    # basic get/add test
    seq = random_sequence()
    graph.insert_sequence(seq)
    graph2 = graph.clone()
    
    assert graph.n_unique() > 0
    assert graph2.n_unique() == 0


@using(ksize=[21, 31, 41], length=[50000, 500000])
@pytest.mark.benchmark(group='dbg-sequence')
def test_pdbg_n_unique(random_sequence, ksize, length, benchmark, partitioned_graph, store):
    from goetia.utils import check_trait
    graph = partitioned_graph

    sequence = random_sequence()
    kmer_set = set((graph.hash(kmer) for kmer in kmers(sequence, ksize)))
    benchmark(graph.insert_sequence, sequence)

    for kmer in kmer_set:
        assert graph.query(kmer)

    if check_trait(libgoetia.storage.is_probabilistic, type(store)):
        assert abs(len(kmer_set) - graph.n_unique()) < length * .001
    else:
        assert len(kmer_set) == graph.n_unique()
        


@using(ksize=[21, 31, 41], length=[50000, 500000])
@pytest.mark.benchmark(group='dbg-sequence')
def test_pdbg_get_counts(random_sequence, ksize, benchmark, partitioned_graph):
    sequence = random_sequence()
    graph = partitioned_graph
    graph.insert_sequence(sequence)

    counts = benchmark(graph.query_sequence, sequence)
    assert all((count > 0 for count in counts))
