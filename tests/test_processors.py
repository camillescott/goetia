# goetia/tests/test_processors.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest

import csv

from .utils import *
from goetia.dbg import dBG
from goetia.parsing import read_fastx

def test_dbg_inserter(graph, datadir, ksize):
    consumer = type(graph).Processor.build(graph, 10000)
    rfile = datadir('random-20-a.fa')

    n_reads = consumer.process(rfile)

    graph2 = graph.shallow_clone()
    for record in read_fastx(rfile):
        for kmer in kmers(record.sequence, ksize):
            assert graph.get(kmer)
            assert not graph2.get(kmer)
            assert graph.get(kmer) != graph2.get(kmer)

    for record in read_fastx(rfile):
        graph2.insert_sequence(record.sequence)
        for kmer in kmers(record.sequence, ksize):
            assert graph.get(kmer) == graph2.get(kmer)


def test_chunked_dbg_inserter(graph, datadir, ksize):
    consumer = type(graph).Processor.build(graph, 10000)
    rfile = datadir('random-20-a.fa')

    for _ in consumer.chunked_process(rfile):
        pass

    graph2 = graph.shallow_clone()
    for record in read_fastx(rfile):
        for kmer in kmers(record.sequence, ksize):
            assert graph.get(kmer)
            assert not graph2.get(kmer)
            assert graph.get(kmer) != graph2.get(kmer)

    for record in read_fastx(rfile):
        graph2.insert_sequence(record.sequence)
        for kmer in kmers(record.sequence, ksize):
            assert graph.get(kmer) == graph2.get(kmer)


@using(ksize=[21,33], length=100)
def test_processor_kmer_timing(graph, ksize, length, random_fasta):
    N = 1000
    n_kmers = N * (length - ksize + 1)
    interval = 10000

    consumer = type(graph).Processor.build(graph, interval)
    sequences, fasta = random_fasta(N)

    prev_time = 0
    for n_seqs, time, n_skipped in consumer.chunked_process(fasta):
        print(prev_time, time, n_seqs)
        if n_seqs != N:
            assert time >= prev_time + interval
        assert time < (prev_time + interval + length)
        prev_time = time
    assert n_seqs == N
    assert time == n_kmers
