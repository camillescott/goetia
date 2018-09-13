# boink/tests/utils.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
import random

from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import reverse_complement as revcomp

from debruijnal_enhance_o_tron.sequence import *
from debruijnal_enhance_o_tron.fixtures.sequence import (using_ksize,
                                                         using_length)

def counting_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('graph_type', 
                                       ['_ByteStorage'],
                                       indirect=['graph_type'],
                                       ids=lambda t: t)(fixture_func)
    return wrapped


def presence_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('graph_type', 
                                       ['_BitStorage', '_SparseppSetStorage'],
                                       indirect=['graph_type'],
                                       ids=lambda t: t)(fixture_func)
    return wrapped


def oxli_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('graph_type', 
                                       ['_BitStorage', '_ByteStorage'],
                                       indirect=['graph_type'],
                                       ids=lambda t: t)(fixture_func)
    return wrapped
