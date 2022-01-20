#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : test_filters.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 13.03.2020

import filecmp
import os

from .utils import run_shell_cmd


def test_solid_filter(datadir, tmpdir):
    with tmpdir.as_cwd():
        rfile = datadir('random-20-a.fa')
        outfile = 'out.fa'

        solid_cmd = ['goetia',
                     'filter',
                     'solid',
                     '-i',
                     '/dev/stdin',
                     '--pairing-mode',
                     'single', 
                     '-K',
                     '31',
                     '-x',
                     '1e7',
                     '-C',
                     '2',
                     '-P',
                     '.9',
                     '-o',
                     outfile]

        ret = run_shell_cmd(['cat', rfile, '|'] + solid_cmd)
        assert ret.returncode == 0
        print(ret.stderr)
        result = open(outfile).read()
        print(f'Result: {result}')
        assert not os.path.getsize(outfile)

        run_shell_cmd(['cat', rfile, rfile, '|'] + solid_cmd)
        print(f'Result: {open(outfile).read()}')
        assert filecmp.cmp(outfile, rfile)

