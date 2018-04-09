# boink/tests/test_processors.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest

import csv

from boink.tests.utils import *

from khmer._oxli.parsing import FastxParser
from boink.cdbg import StreamingCompactor
from boink.processors import FileConsumer, DecisionNodeProcessor


@pytest.mark.parametrize('graph_type', ['BitStorage'], indirect=['graph_type'])
def test_fileconsumer(graph, datadir, ksize):
    consumer = FileConsumer(graph)
    rfile = datadir('random-20-a.fa')

    n_reads, n_kmers = consumer.process(rfile)

    rparser = FastxParser(rfile)
    graph2 = graph.clone()
    for record in rparser:
        for kmer in kmers(record.sequence, ksize):
            assert not graph2.get(kmer)
            assert graph.get(kmer) != graph2.get(kmer)

    for recoord in FastxParser(rfile):
        graph2.add_sequence(record.sequence)
        for kmer in kmers(record.sequence, ksize):
            assert graph.get(kmer) == graph2.get(kmer)


@pytest.mark.parametrize('graph_type', ['BitStorage'], indirect=['graph_type'])
def test_DecisionNodeProcessor(graph, ksize, right_fork, fastx_writer, tmpdir):
    '''TODO Check for false positives
    '''
    sequences = []
    for subgraph_n in range(100):
        (core, branch), _ = right_fork()
        sequences.extend((branch, core))
    fastx_file = fastx_writer(sequences)

    compactor = StreamingCompactor(graph)
    result_file = tmpdir.join('test.csv')
    consumer = DecisionNodeProcessor(compactor, str(result_file))

    n_reads, n_kmers = consumer.process(str(fastx_file))

    data = []
    with open(result_file) as csvfp:
        reader = csv.DictReader(csvfp, skipinitialspace=True)
        for row in reader:
            data.append(row)

            assert int(row['l_degree']) == 1
            assert int(row['r_degree']) == 2

    assert len(data) == 100
