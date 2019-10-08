# boink/tests/test_hashing.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from .utils import *
from boink import libboink
#from boink.hashing import RollingHashShifter, UKHShifter, unikmer_valid


def get_min_unikmer(wmer, uk_map):
    wmer_hash = libboink.hashing.hash_cyclic(wmer, len(wmer))
    kmer_hashes = [libboink.hashing.hash_cyclic(kmer, uk_map.K) for kmer in kmers(wmer, uk_map.K)]
    unikmers = []
    for h in kmer_hashes:
        unikmer = libboink.hashing.UKHS.Unikmer(h)
        if uk_map.query(unikmer):
            unikmers.append(unikmer)
    #print(kmer_hashes)
    #print([str(u) for u in unikmers])
    return wmer_hash, min(unikmers, key=lambda elem: elem.hash)


@pytest.fixture
def unikmer_shifter(request, ksize):
    m = load_unikmer_map(ksize, 7)
    return libboink.hashing.UKHS.LazyShifter(ksize, 7, m), 7, m


def test_rolling_hash():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    hasher = libboink.hashing.RollingHashShifter(K)
    assert hasher.hash(seq) == 13194817695400542713


def test_rolling_hash_seqcursor_eq():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    hasher = libboink.hashing.RollingHashShifter(K)
    hasher.set_cursor(seq)

    assert hasher.hash(seq) == hasher.get()

def test_rolling_setcursor_seq_too_small():
    hasher = libboink.hashing.RollingHashShifter(20)
    with pytest.raises(Exception):
        hasher.set_cursor('AAAAA')

def test_rolling_hash_seq_too_small():
    hasher = libboink.hashing.RollingHashShifter(20)
    with pytest.raises(Exception):
        hasher.hash('AAAAA')


def test_rolling_hash_seq_too_large():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGCAA'
    hasher = libboink.hashing.RollingHashShifter(K)

    assert hasher.hash(seq) == 13194817695400542713


def test_rolling_setcursor_seq_too_large():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGCAA'
    hasher = libboink.hashing.RollingHashShifter(K)

    hasher.set_cursor(seq)
    assert hasher.get() == 13194817695400542713


def test_unikmer_shifter_kmeriterator(ksize, length, random_sequence, unikmer_shifter):
    print()
    shifter, uk_ksize, uk_map = unikmer_shifter
    seq = random_sequence()
    it = libboink.hashing.KmerIterator[type(shifter)](seq, shifter)
    i = 0
    while not it.done():
        kmer = seq[i:i+ksize]
        act_h = it.next()
        #print(i, act_h)
        exp_kmer_hash, exp_ukmer = get_min_unikmer(kmer, uk_map)

        assert act_h.hash == exp_kmer_hash
        assert act_h.unikmer.hash == exp_ukmer.hash
        assert act_h.unikmer.partition == exp_ukmer.partition
        i += 1

        #print('---')


'''
def test_ukhs_unikmer():
    W = 27
    K = 7
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    hasher = UKHShifter(27, 7)
    assert hasher.hash(seq) == 13194817695400542713
    hasher.set_cursor(seq)
    assert hasher.hashvalue == 13194817695400542713
    assert hasher.unikmers() == [(5571541805904823269, 1681)]


@using_length(1000)
def test_ukhs_long_list(linear_path):
    W = 27
    K = 7
    seq = linear_path()

    hasher = UKHShifter(27, 7)
    print(seq[:27])
    U = hasher.find_unikmers(seq)
    print(len(U))
    for u, p in U:
        assert unikmer_valid(p)
'''
