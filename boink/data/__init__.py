#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : __init__.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 06.08.2019

import gzip
import os

from boink.metadata import DATA_DIR

from boink import boink as libboink
from cppyy.gbl import std

def parse_unikmers(W, K):
    import gzip

    valid_W = list(range(20, 210, 10))
    valid_K = list(range(7, 11))
    W = W - (W % 10)

    if not W in valid_W:
        raise ValueError('Invalid UKHS window size.')
    if not K in valid_K:
        raise ValueError('Invalid UKHS K.')

    filename = os.path.join(DATA_DIR,
                            'res_{0}_{1}_4_0.txt.gz'.format(K, W))
    unikmers = std.vector[std.string]()
    with gzip.open(filename, 'rt') as fp:
        for line in fp:
            unikmers.push_back(line.strip())

    return unikmers


def load_unikmer_map(W, K):
    unikmers = parse_unikmers(W, K)
    ukhs = std.make_shared[libboink.hashing.UKHS.Map](W, K, unikmers)
    return ukhs

