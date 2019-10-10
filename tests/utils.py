# boink/tests/utils.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
import collections
import six

import random

from debruijnal_enhance_o_tron.sequence import *
from debruijnal_enhance_o_tron.fixtures.sequence import (using_ksize,
                                                         using_length,
                                                         using_pivot)
#from boink import boink as libboink
from cppyy.gbl import std
import boink
from boink import boink as libboink
from boink.data import load_unikmer_map
from boink.hashing import _types as hashing_types
from boink.utils import check_trait
from boink.storage import _types as storage_types


def storage_t_name(t):
    return t[0].__name__


@pytest.fixture(params=storage_types, ids=storage_t_name)
def storage_type(request):
    _storage_type, params = request.param
    return _storage_type, params


@pytest.fixture(params=hashing_types, ids=lambda t: t.__name__)
def hasher_type(request, ksize):
    _hasher_type = request.param
    if _hasher_type is libboink.hashing.UKHS.LazyShifter:
        return _hasher_type, (ksize, 7, load_unikmer_map(ksize, 7))
    else:
        return _hasher_type, (ksize,)

@pytest.fixture
def hasher(request, hasher_type, ksize):
    _hasher_type, params = hasher_type
    hasher = _hasher_type(*params)
    hasher.set_cursor('A' * ksize)
    return hasher


@pytest.fixture
def store(storage_type):
    _storage_type, params = storage_type
    storage = _storage_type.build(*params)
    return storage


@pytest.fixture()
def graph(store, hasher, ksize):

    _graph_type = libboink.dBG[type(store), type(hasher)]
    
    return _graph_type.build(hasher, store)


def counting_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [(type, args) for type, args in storage_types \
                                        if check_trait(libboink.storage.is_counting, type)],
                                       indirect=['storage_type'],
                                       ids=storage_t_name)(fixture_func)
    return wrapped


def presence_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [(type, args) for type, args in storage_types \
                                        if not check_trait(libboink.storage.is_counting, type)],
                                       indirect=['storage_type'],
                                       ids=storage_t_name)(fixture_func)
    return wrapped


def exact_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [(type, args) for type, args in storage_types \
                                        if not check_trait(libboink.storage.is_probabilistic, type)],
                                       indirect=['storage_type'],
                                       ids=storage_t_name)(fixture_func)
    return wrapped


def is_iterable(arg):
    return (
        isinstance(arg, collections.abc.Iterable) 
        and not isinstance(arg, six.string_types)
    )


def using(**kwargs):
    def wrapped(fixture_func):
        for param, value in kwargs.items():
            if is_iterable(value):
                value = list(value)
                ids = ['{0}={1}'.format(param, v) for v in value]
            else:
                ids = ['{0}={1}'.format(param, value)]
                value = [value]
            
            fixture_func = pytest.mark.parametrize(param,
                                                   value,
                                                   indirect=True,
                                                   ids=ids)(fixture_func)

        return fixture_func

    return wrapped

