import pytest

import khmer
from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp


@pytest.fixture(params=[50000, 500000, 50000000],
                ids=['small', 'medium', 'large'])
def capacity_args(request):
    return optimal_fp(request.param, 0.0001)


@pytest.fixture(params=[khmer.Nodegraph, khmer.Countgraph],
                ids=['(Type=Nodegraph)', '(Type=Countgraph)'])
def graph(request, capacity_args):

    print('Graph Params:', capacity_args)

    return request.param(K, capacity_args.htable_size, 
                         capacity_args.num_htables)


