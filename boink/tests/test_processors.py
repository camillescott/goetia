# boink/tests/test_processors.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from boink.tests.utils import *

from khmer._oxli.parsing import FastxParser
from boink.processors import FileConsumer


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



