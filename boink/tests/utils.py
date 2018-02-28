import pytest
import random

from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import reverse_complement as revcomp
from khmer.tests.graph_structure_fixtures import *


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


@pytest.fixture(params=[khmer.Nodegraph, khmer.Countgraph],
                ids=['(Type=Nodegraph)', '(Type=Countgraph)'])
def graph(request, capacity_args, K):

    print('Graph Params:', capacity_args)

    return request.param(K, capacity_args.htable_size, 
                         capacity_args.num_htables)



