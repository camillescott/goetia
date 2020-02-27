# goetia/tests/utils.py
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.


from distutils import dir_util
import os

from debruijnal_enhance_o_tron.fixtures.subgraphs import *
from debruijnal_enhance_o_tron.fixtures.sequence import *
from debruijnal_enhance_o_tron.fixtures.collectors import *

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


@pytest.fixture
def random_fasta(random_sequence, fastx_writer):

    def get(N):
        sequences = [random_sequence() for i in range(N)]
        return sequences, str(fastx_writer(sequences))
    
    return get