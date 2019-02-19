# boink/tests/test_hashing.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
from boink.tests.utils import *
from boink.hashing import RollingHashShifter, UKHShifter, unikmer_valid


def test_rolling_hash():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    hasher = RollingHashShifter(K)
    assert hasher.hash(seq) == 13194817695400542713


def test_rolling_hash_seqcursor_eq():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGC'

    hasher = RollingHashShifter(K)
    hasher.set_cursor(seq)

    assert hasher.hash(seq) == hasher.hashvalue

def test_rolling_setcursor_seq_too_small():
    hasher = RollingHashShifter(20)
    with pytest.raises(ValueError):
        hasher.set_cursor('AAAAA')

def test_rolling_hash_seq_too_small():
    hasher = RollingHashShifter(20)
    with pytest.raises(ValueError):
        hasher.hash('AAAAA')


def test_rolling_hash_seq_too_large():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGCAA'
    hasher = RollingHashShifter(K)

    assert hasher.hash(seq) == 13194817695400542713


def test_rolling_setcursor_seq_too_large():
    K = 27
    seq = 'TCACCTGTGTTGTGCTACTTGCGGCGCAA'
    hasher = RollingHashShifter(K)

    hasher.set_cursor(seq)
    assert hasher.hashvalue == 13194817695400542713


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
    U = hasher.find_unikmers(seq)
    for u, p in U:
        assert unikmer_valid(p)
    print(U)
    print(len(U))
