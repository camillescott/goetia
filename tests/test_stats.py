import pytest

#from .utils import *
import khmer
from khmer._oxli.hashing import Kmer
from khmer._oxli.parsing import Sequence

from boink.stats import KmerCountFunction, KmerDegreeFunction


def test_KmerCountFunction():
    kmer = Kmer('AAAAA')
    g = khmer.Nodegraph(5, 1e4, 4)
    g.add(str(kmer))
    f = KmerCountFunction(g)
    assert f.evaluate_kmer(kmer) == 1

def test_KmerDegreeFunction():
    kmer = Kmer('ACCTA')
    g = khmer.Nodegraph(5, 1e4, 4)
    g.add(str(kmer))
    f = KmerDegreeFunction(g)
    assert f.evaluate_kmer(kmer) == 0


