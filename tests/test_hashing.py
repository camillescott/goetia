# boink/tests/test_hashing.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from .utils import *
from boink import libboink
from boink.hashing import (FwdLemireShifter, CanLemireShifter, 
                           FwdUnikmerShifter, CanUnikmerShifter,
                           extender_selector_t)

known_kmer = 'TCACCTGTGTTGTGCTACTTGCGGCGC'
known_fwd = 13194817695400542713
known_rc = 4324216031038051805
known_can = known_rc
known_hashes = (known_fwd, known_rc)


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
    return wmer_hash, min(unikmers, key=lambda elem: elem.value)

@using(ksize=27)
def test_extender_hash_base_override(hasher, ksize):
    extender = extender_selector_t[type(hasher)](hasher)

    extender.hash_base(known_kmer)
    assert extender.get_cursor() == known_kmer
    assert extender.get().value in known_hashes


@using(ksize=27)
def test_extender_set_cursor(hasher, ksize):
    extender = extender_selector_t[type(hasher)](hasher)

    extender.set_cursor(known_kmer)
    assert extender.get_cursor() == known_kmer
    assert extender.get().value in known_hashes


@using(ksize=27)
@pytest.mark.parametrize('hasher_type', [FwdLemireShifter, CanLemireShifter], indirect=True)
def test_rolling_hash(hasher, ksize):
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    h = hasher.hash(seq)
    e = 13194817695400542713
    if h.value != e:
        assert e in (h.fw_hash, h.rc_hash)
    else:
        assert h.value == e


@using(ksize=21, length=10000)
def test_hashing_models_hash(hasher, ksize, random_sequence):
    '''test that __hash__ on Hash and Canonical
    yields the value attribute.
    '''

    seq = random_sequence()
    seq += hasher.alphabet.reverse_complement(seq)

    hash_set = set()
    value_set = set()
    for kmer in kmers(seq, ksize):
        h = hasher.hash(kmer)
        # doing hash(h.value) is on purpose! because CPython silently transforms
        # the result of __hash__ into a 62 bit signed range... sigh
        assert hash(h) == hash(h.value), (h, h.value)
        hash_set.add(h)
        value_set.add(h.value)
    
    assert len(hash_set) == len(value_set)
    assert sorted((h.value for h in hash_set)) == sorted((h for h in value_set))


@using(ksize=27)
def test_rolling_nonvolatile_hash(hasher, ksize):
    '''test that hasher.hash([str]) on an instance does not
    change the hasher state.
    '''
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


@using(ksize=27)
def test_hash_base_iter_overload(hasher):
    seq = std.string('TCACCTGTGTTGTGCTACTTGCGGCGC')

    h1 = hasher.hash_base(seq.begin(), seq.end())
    h2 = hasher.hash(seq)

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
        if not hasattr(h, 'fw_hash'):
            _h = h.hash
        else:
            _h = h
        assert e in (_h.fw_hash, _h.rc_hash)
    else:
        assert h.value == e


def test_canonical_rolling_hash(ksize, length, random_sequence):
    seq = random_sequence()

    can_hasher = libboink.hashing.CanLemireShifter(ksize)
    fwd_hasher = libboink.hashing.FwdLemireShifter(ksize)

    for kmer in kmers(seq, ksize):
        rc_kmer = fwd_hasher.alphabet.reverse_complement(kmer)
        fw, rc = fwd_hasher.hash(kmer), fwd_hasher.hash(rc_kmer)
        can = fw if fw < rc else rc

        assert can_hasher.hash(kmer).value == can.value


@pytest.mark.parametrize('hasher_type', [FwdUnikmerShifter, CanUnikmerShifter], indirect=True)
def test_unikmer_hash_base(ksize, length, random_sequence, hasher):
    seq = random_sequence()
    hasher = hasher.build(ksize, 7)

    for i, kmer in enumerate(kmers(seq, ksize)):
        h = hasher.hash_base(kmer)
        exp_kmer_hash, exp_ukmer = get_min_unikmer(kmer, hasher.ukhs_map, type(hasher).base_shifter_type)
        assert h.value == exp_kmer_hash.value, (i, h)
        assert h.minimizer == exp_ukmer, (i, h)


@pytest.mark.parametrize('hasher_type', [FwdUnikmerShifter, CanUnikmerShifter], indirect=True)
def test_unikmer_shift_right(ksize, length, random_sequence, hasher):
    seq = random_sequence()

    h = hasher.hash_base(seq[:ksize])
    exp_kmer_hash, exp_ukmer = get_min_unikmer(seq[:ksize], hasher.ukhs_map, type(hasher).base_shifter_type)
    assert h.value == exp_kmer_hash.value
    assert h.minimizer == exp_ukmer

    for i in range(1, len(seq) - ksize):
        h = hasher.shift_right(seq[i-1], seq[i + ksize - 1])
        exp_hash, exp_uk = get_min_unikmer(seq[i:i+ksize], hasher.ukhs_map, type(hasher).base_shifter_type)

        assert h.value == exp_hash.value
        assert h.minimizer == exp_uk


@pytest.mark.parametrize('hasher_type', [FwdUnikmerShifter, CanUnikmerShifter], indirect=True)
def test_unikmer_shift_left(ksize, length, random_sequence, hasher):
    seq = random_sequence()

    h = hasher.hash_base(seq[-ksize:])
    exp_kmer_hash, exp_ukmer = get_min_unikmer(seq[-ksize:], hasher.ukhs_map, type(hasher).base_shifter_type)
    assert h.value == exp_kmer_hash.value
    assert h.minimizer == exp_ukmer

    for i in range(len(seq) - ksize - 1, -1, -1):
        h = hasher.shift_left(seq[i], seq[i+ksize])
        exp_hash, exp_uk = get_min_unikmer(seq[i:i+ksize], hasher.ukhs_map, type(hasher).base_shifter_type)

        assert h.value == exp_hash.value
        assert h.minimizer == exp_uk


@using(length=30, ksize=27)
def test_shift_right(hasher, ksize, length, random_sequence):
    s = random_sequence()
    
    exp = [hasher.hash(kmer).value for kmer in kmers(s, ksize)]
    fwd_hashes = [hasher.hash_base(s[:ksize]).value]

    for out_base, in_base in zip(s[0:-ksize], s[ksize:]):
        fwd_hashes.append(hasher.shift_right(out_base, in_base).value)

    assert exp == fwd_hashes


@pytest.mark.parametrize('hasher_type', [FwdLemireShifter, CanLemireShifter], indirect=True)
@using(ksize=27)
def test_kmeriterator_owner_init(hasher_type, ksize):
    hasher_type, _ = hasher_type
    it = libboink.hashing.KmerIterator[hasher_type](known_kmer, ksize)
    assert it.first().value in known_hashes


#@pytest.mark.parametrize('hasher_type', [FwdLemireShifter, CanLemireShifter], indirect=True)
#@using(ksize=27)
#def test_kmeriterator_nonowner_init(hasher, ksize, length, random_sequence):
#    it = libboink.hashing.KmerIterator[type(hasher)](known_kmer, hasher.__smartptr__().get())
#    assert it.first().value in known_hashes
#    assert hasher.get().value in known_hashes


@pytest.mark.parametrize('hasher_type', [FwdLemireShifter, CanLemireShifter], indirect=True)
@using(ksize=27)
def test_kmeriterator_proto_init(hasher, ksize):
    it = libboink.hashing.KmerIterator[type(hasher)](known_kmer, hasher)
    assert it.first().value in known_hashes


@pytest.mark.parametrize('hasher_type', [FwdLemireShifter, CanLemireShifter], indirect=True)
def test_kmeriterator(hasher, ksize, length, random_sequence):
    s = random_sequence()

    exp = [hasher.hash(kmer).value for kmer in kmers(s, ksize)]
    
    it = libboink.hashing.KmerIterator[type(hasher)](s, ksize)
    act = []
    while not it.done():
        h = it.next()
        act.append(h.value)
    
    assert act == exp


@pytest.mark.parametrize('hasher_type', [FwdLemireShifter, CanLemireShifter], indirect=True)
def test_kmeriterator_from_proto(hasher, ksize, length, random_sequence):
    s = random_sequence()

    exp = [hasher.hash(kmer).value for kmer in kmers(s, ksize)]
    
    it = libboink.hashing.KmerIterator[type(hasher)](s, hasher)
    act = []
    while not it.done():
        h = it.next()
        act.append(h.value)
    
    assert act == exp



@pytest.mark.parametrize('hasher_type', [FwdLemireShifter, CanLemireShifter], indirect=True)
def test_kmeriterator_hashextender(hasher, ksize, length, random_sequence):
    s = random_sequence()
    extender = extender_selector_t[type(hasher)](hasher)

    exp = [extender.hash(kmer).value for kmer in kmers(s, ksize)]
    
    it = libboink.hashing.KmerIterator[type(extender)](s, ksize)
    act = []
    while not it.done():
        h = it.next()
        act.append(h.value)
    
    assert act == exp

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
