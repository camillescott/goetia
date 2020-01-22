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


@using(ksize=27)
def test_rolling_hash(hasher, ksize):
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    h = hasher.hash(seq)
    e = 13194817695400542713
    if h.value != e:
        assert e in (h.fw_hash, h.rc_hash)
    else:
        assert h.value == e


@using(ksize=27)
def test_rolling_nonvolatile_hash(hasher, ksize):
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    original = hasher.hash_base(seq)
    hasher.hash('A' * ksize)
    assert hasher.get() == original


@using(ksize=27)
def test_hash_base_eq(hasher):
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'
    h1 = hasher.hash(seq)
    h2 = hasher.hash_base(seq)

    assert h1 == h2


@using(ksize=27)
def test_hash_base_get(hasher):
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    h1 = hasher.hash_base(seq)
    h2 = hasher.get()

    assert h1 == h2


@pytest.mark.parametrize('hash_method', ['hash_base', 'hash'])
def test_hash_seq_too_small(hasher, hash_method):
    with pytest.raises(Exception):
        getattr(hasher, hash_method)('AAAAA')


def test_static_hash_seq_too_small(hasher_type):
    with pytest.raises(Exception):
        hasher_type.hash('AAAAA', 7)


@pytest.mark.parametrize('hash_method', ['hash_base', 'hash'])
@using(ksize=27)
def test_seq_too_large(hasher, ksize, hash_method):
    '''A sequence of length > K should have return the
    hash value for sequence[:K]'''

    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGCAA'

    h = getattr(hasher, hash_method)(seq)
    e = 13194817695400542713
    if h.value != e:
        assert e in (h.fw_hash, h.rc_hash)
    else:
        assert h.value == e


def test_canonical_rolling_hash(ksize, length, random_sequence):
    seq = random_sequence()

    can_hasher = libboink.hashing.CanRollingShifter(ksize)
    fwd_hasher = libboink.hashing.FwdRollingShifter(ksize)

    for kmer in kmers(seq, ksize):
        rc_kmer = fwd_hasher.alphabet.reverse_complement(kmer)
        fw, rc = fwd_hasher.hash(kmer), fwd_hasher.hash(rc_kmer)
        can = fw if fw < rc else rc

        assert can_hasher.hash(kmer).value == can.value


'''
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
'''


@using(length=30, ksize=27)
def test_shift_right(hasher, ksize, length, random_sequence):
    s = random_sequence()
    
    exp = [hasher.hash(kmer).value for kmer in kmers(s, ksize)]
    fwd_hashes = [hasher.hash_base(s[:ksize]).value]

    for out_base, in_base in zip(s[0:-ksize], s[ksize:]):
        fwd_hashes.append(hasher.shift_right(out_base, in_base).value)

    assert exp == fwd_hashes


@using(length=30, ksize=27)
def test_shift_right_left_right(hasher, ksize, length, random_sequence):
    s = random_sequence()
    print(s)

    fwd_hashes = [hasher.hash_base(s[:ksize]).value]

    for out_base, in_base in zip(s[0:-ksize], s[ksize:]):
        print(out_base, in_base)
        fwd_hashes.append(hasher.shift_right(out_base, in_base).value)


    bkw_hashes = [hasher.hash_base(s[-ksize:]).value]

    for in_base, out_base in zip(s[:-ksize][::-1], s[ksize:][::-1]):
        print(in_base, out_base)
        bkw_hashes.append(hasher.shift_left(in_base, out_base).value)

    assert fwd_hashes == bkw_hashes[::-1]

    # cursor is now the front k-mer, switch directions and hash fwd again
    # to be sure left->right direction change works
    hashes = [hasher.get().value]
    for out_base, in_base in zip(s[0:-ksize], s[ksize:]):
        hashes.append(hasher.shift_right(out_base, in_base).value)

    assert hashes == fwd_hashes


'''
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
'''
