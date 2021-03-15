#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : test_signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 16.10.2019

import random

from goetia.hashing import Canonical, StrandAware
from goetia.parsing import read_fastx
from goetia.sketches import SourmashSketch, UnikmerSketch, HLLCounter
from goetia.storage import SparseppSetStorage

from .utils import *

import cppyy
import numpy as np
import pandas as pd
import pytest


def test_sourmash_signature(datadir, ksize):
    import sourmash

    rfile = datadir('random-20-a.fa')

    goetia_sig = SourmashSketch.Sketch.build(10000, 31, False, False, False, 42, 0)
    sourmash_sig = sourmash.MinHash(10000, 31)

    processor = SourmashSketch.Processor.build(goetia_sig)
    processor.process(rfile)

    for record in read_fastx(rfile):
        sourmash_sig.add_sequence(record.sequence)

    goetia_mh = goetia_sig.to_sourmash()

    assert goetia_mh.similarity(sourmash_sig) == 1.0


def test_sourmash_scaled(datadir, ksize):
    import sourmash

    rfile = datadir('random-20-a.fa')
    goetia_sig = SourmashSketch.Sketch.build(0, 31, False, False, False, 42, 1000)
    sourmash_sig = sourmash.MinHash(0, 31, scaled=1000)

    processor = SourmashSketch.Processor.build(goetia_sig)
    processor.process(rfile)

    for record in read_fastx(rfile):
        sourmash_sig.add_sequence(record.sequence)

    goetia_mh = goetia_sig.to_sourmash()

    assert goetia_mh.similarity(sourmash_sig) == 1.0


def test_draff_to_numpy(datadir):
    rfile = datadir('random-20-a.fa')

    sketch_t = UnikmerSketch[SparseppSetStorage, StrandAware]
    sketch = sketch_t.Sketch.build(31, 7)
    processor = sketch_t.Processor.build(sketch)
    processor.process(rfile)

    np_sig = sketch.to_numpy()
    py_sig = list(sketch.get_sketch_as_vector())

    assert len(np_sig) == len(py_sig)
    assert list(np_sig) == list(py_sig)

    for np_val, py_val in zip(list(np_sig), list(py_sig)):
        assert np_val == py_val


@pytest.mark.xfail()
def test_hllcounter():
    ints = set(np.random.randint(0, 100000, 100000))
    e = 0.01
    hll = HLLCounter(e)

    for i in ints:
        hll.insert(cppyy.gbl.uint64_t(i))

    est = hll.estimate_cardinality()
    act = len(ints)
    
    margin = (2 * e) * act

    print(act)
    print(est)
    print(hll.get_ncounters())
    print(hll.get_error_rate())
    print(hll.get_alpha())
    print(hll.get_erate())

    assert (act - margin) < est < (act + margin)


def test_sourmash_stream(tmpdir, datadir):
    with tmpdir.as_cwd():
        rfile = datadir('sacPom.pombase.fa.gz')

        goetia_sig = 'goetia.sig'
        sourmash_sig = 'sourmash.sig'
        comp_file = 'comp.csv'

        goetia_cmd = ['goetia', 'sketch', 'sourmash',
                      '--interval', '50000', 
                      '-K', '31', '-N', '100000', '--save-sig', goetia_sig,
                      '-i', rfile]
        sourmash_cmd = ['sourmash', 'compute', '-k', '31', '-n', '100000',
                        '-o', sourmash_sig, rfile]
        compare_cmd = ['sourmash', 'compare', goetia_sig, sourmash_sig, '--csv', comp_file]

        run_shell_cmd(goetia_cmd)
        run_shell_cmd(sourmash_cmd)
        run_shell_cmd(compare_cmd)

        comp = pd.read_csv(comp_file)
        assert comp.all().all()

