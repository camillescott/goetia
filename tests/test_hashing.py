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
    positions = []
    for i, h in enumerate(kmer_hashes):
        unikmer = libboink.hashing.UKHS.Unikmer(h)
        if uk_map.query(unikmer):
            unikmers.append(unikmer)
            positions.append(i)
    print([f'{i}: {u}' for i, u in zip(positions, unikmers)])
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


def test_unikmer_shifter_shift_left(ksize, length, random_sequence, unikmer_shifter):
    shifter, uk_ksize, uk_map = unikmer_shifter
    seq = random_sequence()

    hashes = [shifter.set_cursor(seq[-ksize:])]
    exp_kmer_hash, exp_ukmer = get_min_unikmer(seq[-ksize:], uk_map)
    assert hashes[0].hash == exp_kmer_hash
    assert hashes[0].unikmer == exp_ukmer

    print(seq[-ksize:])
    for i in range(len(seq) - ksize - 1, -1, -1):
        print(seq[i:i+ksize])
        h = shifter.shift_left(seq[i])
        exp_hash, exp_uk = get_min_unikmer(seq[i:i+ksize], uk_map)

        assert h.hash == exp_hash
        assert h.unikmer == exp_uk


def test_update_left_right(hasher, ksize, length, random_sequence):
    s = random_sequence()
    fwd_hashes = [hasher.set_cursor(s[:ksize])]
    for base in s[ksize:]:
        fwd_hashes.append(hasher.shift_right(base))

    bkw_hashes = [hasher.set_cursor(s[-ksize:])]
    for base in s[:-ksize][::-1]:
        bkw_hashes.append(hasher.shift_left(base))

    assert fwd_hashes == bkw_hashes[::-1]

    # cursor is now the front k-mer, switch directions and hash fwd again
    # to be sure left->right direction change works
    hashes = [hasher.get()]
    for base in s[ksize:]:
        hashes.append(hasher.shift_right(base))

    assert hashes == fwd_hashes


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
