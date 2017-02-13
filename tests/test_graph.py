import pytest

from boink.graph import ProbabilisticGraph
from boink.stats import PSequenceFunction, PKmerFunction

from khmer._oxli.hashing import Kmer
from khmer import Nodegraph

def test_basic_kmer_insert():
    '''The default functions always return P=1.0'''

    g = Nodegraph(5, 1e3, 4)
    pg = ProbabilisticGraph(g, 
                            kmer_func=PKmerFunction(),
                            sequence_func=PSequenceFunction())
    kmer = Kmer('AAAAA')
    pg.insert(kmer)
    assert g.get(str(kmer))
