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

import cppyy.gbl
from cppyy.gbl import std
import boink
from boink import boink as libboink
from boink.hashing import types as hashing_types
from boink.utils import check_trait, pretty_repr
from boink.storage import _types as storage_types

from boink.hashing import (FwdRollingShifter, CanRollingShifter,
                           FwdUnikmerShifter, CanUnikmerShifter,
                           UKHS)


@pytest.fixture(params=storage_types, ids=lambda t: pretty_repr(t[0]))
def storage_type(request):
    _storage_type, params = request.param
    return _storage_type, params


@pytest.fixture(params=[FwdRollingShifter, CanRollingShifter,
                        FwdUnikmerShifter, CanUnikmerShifter],
                ids=lambda t: pretty_repr(t))
def hasher_type(request, ksize):
    _hasher_type = request.param
    if 'Unikmer' in _hasher_type.__name__:
        return _hasher_type, (ksize, 7)
    else:
        return _hasher_type, (ksize,)

@pytest.fixture
def hasher(request, hasher_type, ksize):
    _hasher_type, params = hasher_type
    if hasattr(_hasher_type, 'build'):
        hasher = _hasher_type.build(*params)
    else:
        hasher = _hasher_type(*params)
    #hasher.set_cursor('A' * ksize)
    return hasher


@pytest.fixture
def store(storage_type):
    _storage_type, params = storage_type
    storage = _storage_type.build(*params)
    return storage


@pytest.fixture
def partitioned_graph(store, ksize):
    ukhs = UKHS[FwdRollingShifter].load(ksize, 7)
    pstore = std.make_shared[libboink.storage.PartitionedStorage[type(store)]](ukhs.n_hashes(),
                                                                               store)
    graph = std.make_shared[libboink.PdBG[type(store), FwdUnikmerShifter]](ksize, 7,
                                                                           ukhs,
                                                                           pstore)
    return graph


@pytest.fixture()
def graph(store, hasher, ksize):

    _graph_type = libboink.dBG[type(store), type(hasher)]
    
    return _graph_type.build(store, hasher)


def counting_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [(type, args) for type, args in storage_types \
                                        if check_trait(libboink.storage.is_counting, type)],
                                       indirect=['storage_type'],
                                       ids=lambda t: pretty_repr(t[0]))(fixture_func)
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
                                       ids=lambda t: pretty_repr(t[0]))(fixture_func)
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
                                       ids=lambda t: pretty_repr(t[0]))(fixture_func)
    return wrapped


def is_iterable(arg):
    return (
        isinstance(arg, collections.abc.Iterable) 
        and not isinstance(arg, six.string_types)
    )


def using(**kwargs):

    def pretty(val):
        if 'meta' in type(val).__name__:
            return pretty_repr(val)
        else:
            return str(val)

    def wrapped(fixture_func):
        for param, value in kwargs.items():
            if is_iterable(value):
                value = list(value)
                ids = ['{0}={1}'.format(param, pretty(v)) for v in value]
            else:
                ids = ['{0}={1}'.format(param, pretty(value))]
                value = [value]
            
            fixture_func = pytest.mark.parametrize(param,
                                                   value,
                                                   indirect=True,
                                                   ids=ids)(fixture_func)

        return fixture_func

    return wrapped

