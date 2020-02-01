#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : test_ukhs.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 22.01.2020


import pytest
from .utils import *
from boink import libboink
from boink.hashing import (FwdLemireShifter, CanLemireShifter, UKHS)


@using(ksize=[21,31,41])
def test_ukhs_query(ksize):
    utype = UKHS[FwdLemireShifter]
    ukhs = utype.load(ksize, 7)
    unikmers = utype.parse_unikmers(ksize, 7)

    for i, kmer in enumerate(unikmers):
        uhash = FwdLemireShifter.hash(kmer, 7)
        result = ukhs.query(uhash).value()

        assert uhash == result.value
        assert result.partition == i
