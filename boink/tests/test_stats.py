# boink/tests/test_partitioning.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest

#from .utils import *
import khmer
from khmer._oxli.hashing import Kmer
from khmer._oxli.sequence import Sequence

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


