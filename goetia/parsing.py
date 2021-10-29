#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : parsing.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 27.02.2020

from goetia.utils import find_common_basename, remove_fx_suffix, grouper

from goetia import libgoetia
from goetia.alphabets import DNA_SIMPLE

import os

PAIRING_MODES     = ('split', 'interleaved', 'single')

FastxParser       = libgoetia.parsing.FastxParser
SplitPairedReader = libgoetia.parsing.SplitPairedReader


def read_fastx(filename, alphabet=DNA_SIMPLE, strict=False, min_length=0):
    reader = FastxParser[alphabet].build(filename, strict, min_length)

    for record in reader:
        yield record
    

def get_fastx_args(parser, default_mode='single'):
    """Common pairing mode argument."""
    group = parser.add_argument_group('fastx arguments')
    group.add_argument('--pairing-mode', default=default_mode,
                        choices=PAIRING_MODES,
                        help='How to interpret read pairing. With `single`, '\
                             'reads will be parsed as singletons, regardless'\
                             ' of pairing or file order. With `interleaved`,'\
                             ' each file will be assumed to be interleaved '\
                             'and paired, with singletons allowed to be mixed'\
                             ' in. With `split`, it will be assumed that each'\
                             ' group of two files in the input list are '\
                             'as (LEFT, RIGHT), ...')
    group.add_argument('--names', nargs='+',
                       help='Names to associate with samples. If pairing mode '
                            'is "single" or "interleaved," this should have length equal to the '
                            'number of samples; if "split," it should have an entry for each pair.')
    return group


def iter_fastx_inputs(inputs, pairing_mode, names=None, merge=False):
    if pairing_mode == 'split':
        _samples  = list(grouper(2, inputs))
        _names = []
        for n, (l, r) in enumerate(_samples):
            base = find_common_basename(remove_fx_suffix(l),
                                        remove_fx_suffix(r))
            _names.append(base or f'Sample {n+1}')
    else:
        _samples  = [(s, ) for s in inputs]
        _names = [os.path.basename(remove_fx_suffix(s)) for s, in _samples]

    if merge:
        _names = [merge] * len(_samples)

    if names:
        if merge:
            pass
        else:
            if len(names) != len(_samples):
                raise RuntimeError('Number of names must match number of samples.')
            _names = names

    yield from zip(_samples, _names)
