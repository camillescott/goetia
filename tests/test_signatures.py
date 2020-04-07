#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : test_signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.10.2019

import pytest

from goetia.signatures import SourmashSketch

from .utils import *
import screed

def test_sourmash_signature(datadir, ksize):
    import sourmash

    rfile = datadir('random-20-a.fa')

    goetia_sig = SourmashSketch.Signature.build(10000, 31, False, False, False, 42, 0)
    sourmash_sig = sourmash.MinHash(10000, 31)

    processor = SourmashSketch.Processor.build(goetia_sig)
    processor.process(rfile)

    for record in screed.open(rfile):
        sourmash_sig.add_sequence(record.sequence)

    goetia_mh = goetia_sig.to_sourmash()

    assert goetia_mh.similarity(sourmash_sig) == 1.0


def test_sourmash_scaled(datadir, ksize):
    import sourmash

    rfile = datadir('random-20-a.fa')
    goetia_sig = SourmashSketch.Signature.build(0, 31, False, False, False, 42, 1000)
    sourmash_sig = sourmash.MinHash(0, 31, scaled=1000)

    processor = SourmashSketch.Processor.build(goetia_sig)
    processor.process(rfile)

    for record in screed.open(rfile):
        sourmash_sig.add_sequence(record.sequence)

    goetia_mh = goetia_sig.to_sourmash()

    assert goetia_mh.similarity(sourmash_sig) == 1.0


