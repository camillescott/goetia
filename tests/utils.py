# boink/tests/utils.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
import random

from debruijnal_enhance_o_tron.sequence import *
from debruijnal_enhance_o_tron.fixtures.sequence import (using_ksize,
                                                         using_length,
                                                         using_pivot)
#from boink import boink as libboink
from cppyy.gbl import std
from boink import boink as libboink
from boink import types


def storage_t_name(t):
    return t[0].__name__


@pytest.fixture(params=types.storage_types, ids=storage_t_name)
def storage_type(request):
    _storage_type, params = request.param
    return _storage_type, params


@pytest.fixture()
def graph(storage_type, ksize):
    _storage_type, params = storage_type
    _graph_type = libboink.dBG[_storage_type, libboink.hashing.RollingHashShifter]
    
    return std.make_shared[_graph_type](ksize, *params)



# (boink.storage.QFStorage, (20,)),
def counting_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [(libboink.storage.ByteStorage, (1e6, 4)),
                                        (libboink.storage.NibbleStorage, (1e6, 4)),
                                        (libboink.storage.QFStorage, (20,))],
                                       indirect=['storage_type'],
                                       ids=storage_t_name)(fixture_func)
    return wrapped


def presence_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [(libboink.storage.BitStorage, (1e6, 4)),
                                        (libboink.storage.SparseppSetStorage, ())],
                                       indirect=['storage_type'],
                                       ids=storage_t_name)(fixture_func)
    return wrapped


def exact_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [(libboink.storage.SparseppSetStorage, ())],
                                       indirect=['storage_type'],
                                       ids=storage_t_name)(fixture_func)
    return wrapped
