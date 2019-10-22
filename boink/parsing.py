from boink.utils import find_common_basename, remove_fx_suffix, grouper

from boink import libboink

import os

PAIRING_MODES     = ('split', 'interleaved', 'single')

SequenceReader    = libboink.parsing.SequenceReader
FastxParser       = libboink.parsing.FastxParser

FastxReader       = SequenceReader[FastxParser]
SplitPairedReader = libboink.parsing.SplitPairedReader[FastxParser]


def get_pairing_args(parser):
    """Common pairing mode argument."""
    group = parser.add_argument_group('input/output')
    group.add_argument('--pairing-mode', default='interleaved',
                        choices=PAIRING_MODES,
                        help='How to interpret read pairing. With `single`, '\
                             'reads will be parsed as singletons, regardless'\
                             ' of pairing or file order. With `interleaved`,'\
                             ' each file will be assumed to be interleaved '\
                             'and paired, with singletons allowed to be mixed'\
                             ' in. With `split`, it will be assumed that each'\
                             ' group of two files in the input list are '\
                             'as (LEFT, RIGHT), ...')
    return group


def iter_fastx_inputs(inputs, pairing_mode, prefix=None, merge=False):
    if pairing_mode == 'split':
        _samples  = list(grouper(2, inputs))
        _prefixes = [find_common_basename(remove_fx_suffix(l),
                                          remove_fx_suffix(r)) for l, r in _samples]
    else:
        _samples  = [(s, ) for s in inputs]
        _prefixes = [os.path.basename(remove_fx_suffix(s)) for s, in _samples]

    if merge:
        _prefixes = [merge] * len(_samples)

    if prefix:
        if merge:
            pass
        else:
            if len(prefix) != len(_samples):
                raise RuntimeError('Number of names must match number of samples.')
            _prefixes = prefix

    yield from zip(_samples, _prefixes)
