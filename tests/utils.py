# goetia/tests/utils.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest
import collections
import six

import os
import subprocess

from debruijnal_enhance_o_tron.sequence import *

import cppyy.gbl
from cppyy.gbl import std
from goetia import goetia as libgoetia
from goetia.hashing import types as hashing_types
from goetia.utils import pretty_repr
from goetia.storage import types as storage_types

from goetia.hashing import (FwdLemireShifter, CanLemireShifter,
                           FwdUnikmerShifter, CanUnikmerShifter,
                           UKHS)


# don't use HLLStorage for the dBG tests, it isn't meant to work
storage_types = [t for t in storage_types if 'HLLStorage' not in t.__name__]


@pytest.fixture(params=storage_types, ids=lambda t: pretty_repr(t))
def storage_type(request):
    return request.param


@pytest.fixture(params=[FwdLemireShifter, CanLemireShifter,
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
    storage = storage_type.build()
    return storage


@pytest.fixture
def partitioned_graph(storage_type, ksize):

    if storage_type in [libgoetia.SparseppSetStorage,
                        libgoetia.BTreeStorage,
                        libgoetia.PHMapStorage]:
        params = tuple()
    elif storage_type is libgoetia.QFStorage:
        params = (10, )
    else:
        params = (1000, 4)

    store = storage_type.build(*params)
    ukhs = UKHS[FwdLemireShifter].load(ksize, 7)
    pstore = std.make_shared[libgoetia.PartitionedStorage[type(store)]](ukhs.n_hashes(),
                                                                               store)
    graph = std.make_shared[libgoetia.PdBG[type(store), FwdUnikmerShifter]](ksize, 7,
                                                                           ukhs,
                                                                           pstore)
    return graph


@pytest.fixture()
def graph(store, hasher, ksize):

    _graph_type = libgoetia.dBG[type(store), type(hasher)]
    
    return _graph_type.build(store, hasher)


def counting_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [type for type in storage_types if type.is_counting],
                                       indirect=['storage_type'],
                                       ids=lambda t: pretty_repr(t))(fixture_func)
    return wrapped


def presence_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [type for type in storage_types if not type.is_counting],
                                       indirect=['storage_type'],
                                       ids=lambda t: pretty_repr(t))(fixture_func)
    return wrapped


def exact_backends(*args):
    '''
    Convenience wrapper to reduce verbosity of indirect parametrization
    '''
    def wrapped(fixture_func):
        return pytest.mark.parametrize('storage_type', 
                                       [type for type in storage_types if not type.is_probabilistic],
                                       indirect=['storage_type'],
                                       ids=lambda t: pretty_repr(t))(fixture_func)
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


def run_shell_cmd(cmd, in_directory=None):
    cwd = os.getcwd()
    if in_directory:
        os.chdir(in_directory)

    print('running: ', ' '.join(cmd))
    try:
        p = subprocess.run(' '.join(cmd), shell=True, check=True,
                           stderr=subprocess.PIPE)
        return p
    finally:
        os.chdir(cwd)
