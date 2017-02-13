import pytest

#from .utils import *
import khmer
from khmer._oxli.hashing import Kmer

from boink.graph import ProbabilisticGraph
from boink.stats import PSequenceFunction, PKmerFunction, PFunction_shim

def test_PKmerFunc():
    kmer = Kmer('AAAAA')
    g = khmer.Nodegraph(5, 1e4, 4)
    f = PKmerFunction(g)
    assert f.evaluate(kmer) == 1.0

