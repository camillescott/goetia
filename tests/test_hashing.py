# boink/tests/test_hashing.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from .utils import *
from boink import libboink
from boink.hashing import FwdRollingShifter, CanRollingShifter
#from boink.hashing import RollingHashShifter, UKHShifter, unikmer_valid


def get_min_unikmer(wmer, uk_map, shifter_type):
    wmer_hash = shifter_type.hash(wmer, len(wmer))
    kmer_hashes = [shifter_type.hash(kmer, uk_map.K) for kmer in kmers(wmer, uk_map.K)]
    unikmers = []
    positions = []
    for i, h in enumerate(kmer_hashes):
        unikmer = uk_map.query(h)
        if unikmer:
            unikmers.append(unikmer.value())
            positions.append(i)
    print([f'{i}: {u}' for i, u in zip(positions, unikmers)])
    #print([str(u) for u in unikmers])
    return wmer_hash, min(unikmers, key=lambda elem: elem.hash)


def test_fwd_rolling_hash():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    hasher = FwdRollingShifter(K)
    assert hasher.hash(seq).value == 13194817695400542713


@using(ksize=27)
def test_fwd_hash_base_eq(hasher):
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'
    h1 = hasher.hash(seq)
    h2 = hasher.hash_base(seq)

    print(h1, h2)
    assert h1 == h2


@using(ksize=27)
def test_fwd_hash_base_get(hasher):
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    h1 = hasher.hash_base(seq)
    h2 = hasher.get()

    print(h1, h2)

    assert h1 == h2


def test_hash_base_seq_too_small(hasher):
    with pytest.raises(Exception):
        hasher.hash_base('AAAAA')


def test_static_hash_seq_too_small():
    hasher = libboink.hashing.RollingHashShifter(20)
    with pytest.raises(Exception):
        hasher.hash('AAAAA')


def test_fw_rolling_hash_seq_too_large():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGCAA'
    hasher = FwdRollingShifter(K)

    assert hasher.hash(seq).value == 13194817695400542713


def test_rolling_setcursor_seq_too_large():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGCAA'
    hasher = libboink.hashing.RollingHashShifter(K)

    hasher.set_cursor(seq)
    assert hasher.get().value == 13194817695400542713


def test_canonical_rolling_hash(ksize, length, random_sequence):
    seq = random_sequence()

    can_hasher = libboink.hashing.CanRollingShifter(ksize)
    fwd_hasher = libboink.hashing.FwdRollingShifter(ksize)

    for kmer in kmers(seq, ksize):
        rc_kmer = fwd_hasher.alphabet.reverse_complement(kmer)
        fw, rc = fwd_hasher.hash(kmer), fwd_hasher.hash(rc_kmer)
        can = fw if fw < rc else rc

        assert can_hasher.hash(kmer) == can


def test_unikmer_shifter_shift_left(ksize, length, random_sequence, unikmer_shifter):
    shifter, uk_ksize, uk_map = unikmer_shifter
    seq = random_sequence()

    hashes = [shifter.set_cursor(seq[-ksize:])]
    exp_kmer_hash, exp_ukmer = get_min_unikmer(seq[-ksize:], uk_map)
    assert hashes[0].hash == exp_kmer_hash
    assert hashes[0].unikmer == exp_ukmer

    for i in range(len(seq) - ksize - 1, -1, -1):
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

