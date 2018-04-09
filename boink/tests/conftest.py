# boink/tests/utils.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import pytest

from distutils import dir_util
import os

from debruijnal_enhance_o_tron.fixtures.subgraphs import *
from debruijnal_enhance_o_tron.fixtures.sequence import *

from boink.dbg import get_dbg_type


@pytest.fixture(params=['BitStorage', 'ByteStorage', 'NibbleStorage'])
def graph_type(request, ksize):

    _graph_type = get_dbg_type(storage=request.param)

    class BoinkAdapter(_graph_type):
        '''Basic adapter for boink dBG.
        '''

        def __init__(self, ksize, starting_size, n_tables):
            self.ksize = ksize
            self.size = starting_size
            self.n_tables = n_tables
            self.storage = request.param
            pass

        def __new__(cls, ksize, starting_size, n_tables, *args, **kwargs):
            return super().__new__(cls, ksize, starting_size, n_tables,
                                   *args, **kwargs)

        def add(self, item):
            if not isinstance(item, int) and len(item) < self.ksize:
                raise ValueError()
            elif isinstance(item, int) or len(item) == self.ksize:
                return super().add(item)
            else:
                return self.add_sequence(item)

        def __repr__(self):
            return '<BoinkAdapter[{0}] K={1} X={2} N={3}>'.format(self.storage,
                                                                  self.ksize,
                                                                  self.size,
                                                                  self.n_tables)


    return _graph_type, BoinkAdapter


@pytest.fixture
def graph(ksize, graph_type):
    _graph_type, Adapter = graph_type
    return Adapter(ksize, 10000000, 4)


@pytest.fixture(params=[50000, 500000, 50000000],
                ids=['small', 'medium', 'large'])
def capacity_args(request):
    return optimal_fp(request.param, 0.0001)


@pytest.fixture
def fastx_writer(tmpdir):

    def write(sequences):
        filepath = tmpdir.join('tmp.fasta')
        with open(filepath, 'w') as fp:
            for n, sequence in enumerate(sequences):
                fp.write('>{0}\n{1}\n'.format(n, sequence))
        return filepath

    return write


@pytest.fixture
def datadir(tmpdir, request):
    '''
    Fixture responsible for locating the test data directory and copying it
    into a temporary directory.
    '''
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    data_dir = os.path.join(test_dir, 'test-data') 
    dir_util.copy_tree(data_dir, str(tmpdir))

    def getter(filename, as_str=True):
        filepath = tmpdir.join(filename)
        if as_str:
            return str(filepath)
        return filepath

    return getter
