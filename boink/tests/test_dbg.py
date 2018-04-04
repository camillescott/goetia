# boink/tests/test_dbg.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from boink.tests.utils import *


@using_ksize([21, 51, 101])
def test_presence(graph, ksize, random_sequence):
    # basic get/add test
    for kmer in kmers(random_sequence(), ksize):

        hashval = graph.hash(kmer)

        assert graph.get(kmer) == 0
        assert graph.get(hashval) == 0

        graph.add(kmer)
        assert graph.get(kmer) == 1
        assert graph.get(hashval) == 1

        graph.add(kmer)
        # Node* types can only tell presence/absence
        if 'Bit' in graph.storage:
            assert graph.get(kmer) == 1
            assert graph.get(hashval) == 1
        else:
            assert graph.get(kmer) == 2
            assert graph.get(hashval) == 2


@using_ksize([21,151])
def test_n_occupied(graph, ksize):
    # basic get/add test
    kmer = 'G' * ksize

    assert graph.n_occupied == 0
    assert graph.n_unique == 0

    graph.add(kmer)
    assert graph.n_occupied == 1
    assert graph.n_unique == 1

    graph.add(kmer)
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


# TODO add test for k-mer too short

'''

def test_hashsizes(graph):
    # hashsizes method.
    graph = dbg_type(5)
    assert (graph.hashsizes() == PRIMES_1m or
            # CQF allocates some extra slots beyond what you request
            # exactly how many extra is an implementation detail
            graph.hashsizes()[0] >= QF_SIZE)
'''


@using_ksize(5)
def test_add_hashval(graph, ksize):
    # test add(hashval)
    x = graph.hash("ATGGC")
    y = graph.add(x)
    assert y

    z = graph.get(x)
    assert z == 1


@using_ksize(5)
def test_add_dna_kmer(graph, ksize):
    # test add(dna)
    x = graph.add("ATGGC")
    assert x

    z = graph.get("ATGGC")
    assert z == 1


@using_ksize(5)
def test_get_hashval(graph, ksize):
    # test get(hashval)
    hashval = graph.hash("ATGGC")
    graph.add(hashval)

    z = graph.get(hashval)
    assert z == 1


@using_ksize(5)
def test_get_hashval_rc(graph, ksize):
    # fw and rc should NOT be the same on this table
    hashval = graph.hash("ATGC")
    rc = graph.hash("GCAT")

    assert hashval != rc


@using_ksize(5)
def test_get_dna_kmer(graph, ksize):
    # test get(dna)
    hashval = graph.hash("ATGGC")
    graph.add(hashval)

    z = graph.get("ATGGC")
    assert z == 1


@using_ksize(5)
def test_get_bad_dna_kmer(graph, ksize):
    # test get(dna) with bad dna; should fail
    with pytest.raises(ValueError):
        graph.get("ATYGC")


@using_ksize(5)
def test_add_sequence(graph, ksize):
    x = "ATGCCGATGCA"
    num_kmers = sum(graph.add_sequence(x))
    assert num_kmers == len(x) - ksize + 1   # num k-mers consumed

    for start in range(len(x) - 6 + 1):
        assert graph.get(x[start:start + 6]) == 1


@using_ksize(5)
def test_add_sequence_bad_dna(graph):
    # while we don't specifically handle bad DNA, we should at least be
    # consistent...
    x = "ATGCCGNTGCA"
    with pytest.raises(ValueError):
        num_kmers = graph.add_sequence(x)


@using_ksize(10)
def test_consume_short(graph):
    # raise error on too short when consume is run
    x = "ATGCA"
    with pytest.raises(ValueError):
        graph.add_sequence(x)


@using_ksize(6)
def test_get_kmer_counts(graph):
    graph.add_sequence("AAAAAA")
    counts = graph.get_counts("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] == 1

    graph.add_sequence("AAAAAA")
    counts = graph.get_counts("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] >= 1

    graph.add_sequence("AAAAAT")
    counts = graph.get_counts("AAAAAAT")
    print(counts)
    assert len(counts) == 2
    assert counts[0] >= 1
    assert counts[1] == 1


@using_ksize(5)
def test_consume_alias(graph_type, ksize):
    g1 = graph_type(ksize, 1000, 4)
    g2 = graph_type(ksize, 1000, 4)
    assert g1.add_sequence('A' * (ksize + 1)) == \
           g2.consume('A' * (ksize + 1))


@using_ksize(6)
def test_get_kmer_hashes(graph):
    hashes = list(graph.hashes("ACGTGCGT"))
    print(hashes)
    assert len(hashes) == 3
    assert hashes[0] == graph.hash("ACGTGC")
    assert hashes[1] == graph.hash("CGTGCG")
    assert hashes[2] == graph.hash("GTGCGT")


'''
def test_get_min_count(graph):
    graph = dbg_type(6)

    # master string, 3 k-mers
    x = "ACGTGCGT"

    graph.add("ACGTGC")  # 3
    graph.add("ACGTGC")
    graph.add("ACGTGC")

    graph.add("CGTGCG")  # 1

    graph.add("GTGCGT")  # 2
    graph.add("GTGCGT")

    counts = graph.get_kmer_counts(x)
    assert graph.get_min_count(x) == min(counts)
    assert graph.get_max_count(x) == max(counts)
    med, _, _ = graph.get_median_count(x)
    assert med == list(sorted(counts))[len(counts) // 2]


def test_trim_on_abundance(graph):
    graph = dbg_type(6)

    x = "ATGGCAGTAGCAGTGAGC"
    graph.consume(x[:10])

    (y, pos) = graph.trim_on_abundance(x, 1)
    assert pos == 10
    assert x[:pos] == y


def test_trim_below_abundance(graph):
    graph = dbg_type(6)

    x = "ATGGCAGTAGCAGTGAGC"
    x_rc = screed.rc(x)
    graph.consume(x_rc[:10])

    print(len(x))

    (y, pos) = graph.trim_below_abundance(x, 0)
    assert pos == len(x) - graph.ksize() + 1
    assert x[:pos] == y


DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"


def test_consume_seqfile_reads_parser(graph):
    graph = dbg_type(5)
    rparser = ReadParser(utils.get_test_data('test-fastq-reads.fq'))

    graph.consume_seqfile(rparser)

    graph2 = graph(5)
    for record in screed.open(utils.get_test_data('test-fastq-reads.fq')):
        graph2.consume(record.sequence)

    assert graph.get('CCGGC') == graph2.get('CCGGC')


def test_consume_seqfile(dbg_type):
    graph = graph(5)
    graph.consume_seqfile(utils.get_test_data('test-fastq-reads.fq'))

    graph2 = dbg_type(5)
    for record in screed.open(utils.get_test_data('test-fastq-reads.fq')):
        graph2.consume(record.sequence)

    assert graph.get('CCGGC') == graph2.get('CCGGC')


def test_save_load(Tabletype):
    graph = Tabletype(5)
    graphype = type(graph)
    savefile = utils.get_temp_filename('tablesave.out')

    # test add(dna)
    x = graph.add("ATGGC")
    z = graph.get("ATGGC")
    assert z == 1

    graph.save(savefile)

    # should we provide a single load function here? yes, probably. @CTB
    loaded = ttype.load(savefile)

    z = loaded.get('ATGGC')
    assert z == 1


def test_get_bigcount(Tabletype):
    # get_bigcount should return false by default
    tt = Tabletype(12)

    assert not tt.get_use_bigcount()


def test_set_bigcount(Tabletype):
    supports_bigcount = [Countgraph, Counttable, CyclicCounttable]
    tt = Tabletype(12)

    if type(tt) in supports_bigcount:
        tt.set_use_bigcount(True)

        for i in range(300):
            tt.add('G' * 12)
        assert tt.get('G' * 12) == 300

    else:
        with pytest.raises(ValueError):
            tt.set_use_bigcount(True)


def test_abund_dist_A(graph):
    A_filename = utils.get_test_data('all-A.fa')

    graph = dbg_type(4)
    tracking = Nodegraph(4, 1, 1, primes=PRIMES_1m)

    graph.consume_seqfile(A_filename)
    dist = graph.abundance_distribution(A_filename, tracking)

    print(dist[:10])
    assert sum(dist) == 1
    assert dist[0] == 0


def test_abund_dist_A_readparser(graph):
    A_filename = utils.get_test_data('all-A.fa')
    rparser = ReadParser(A_filename)

    graph = dbg_type(4)
    tracking = Nodegraph(4, 1, 1, primes=PRIMES_1m)

    graph.consume_seqfile(A_filename)
    dist = graph.abundance_distribution(rparser, tracking)

    print(dist[:10])
    assert sum(dist) == 1
    assert dist[0] == 0
'''
