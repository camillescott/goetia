import pytest
import random

import khmer
from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import reverse_complement as revcomp


@pytest.fixture(params=[21,25,27,31],
                scope='module')
def K(request):
    return request.param


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


def mutate_base(base):
    if base in 'AT':
        return random.choice('GC')
    elif base in 'GC':
        return random.choice('AT')
    else:
        assert False, 'bad base'


def mutate_sequence(sequence, N=1):
    sequence = list(sequence)
    positions = random.sample(range(len(sequence)), N)

    for i in positions:
        sequence[i] = mutate_base(sequence[i])

    return ''.join(sequence)


def mutate_position(sequence, pos):
    sequence = list(sequence)
    sequence[pos] = mutate_base(sequence[pos])
    return ''.join(sequence)


def test_mutate_sequence():
    for _ in range(100):
        assert 'A' not in mutate_sequence('A' * 10, 10)
        assert 'T' not in mutate_sequence('T' * 10, 10)
        assert 'C' not in mutate_sequence('C' * 10, 10)
        assert 'G' not in mutate_sequence('G' * 10, 10)


def test_mutate_position():
    assert mutate_position('AAAA', 2) in ['AACA', 'AAGA']
    assert mutate_position('TTTT', 2) in ['TTCT', 'TTGT']
    assert mutate_position('CCCC', 2) in ['CCAC', 'CCTC']
    assert mutate_position('GGGG', 2) in ['GGAG', 'GGTG']


def get_random_sequence(length, K, exclude=None):
    '''Generate a random (non-looping) nucleotide sequence.

    To be non-overlapping, the sequence should not include any repeated
    length K-1 k-mers.

    Args:
        exclude (str): If not None, add the k-mers from this sequence to the
        seen set.

    Returns:
        str: A random non-looping sequence.
    '''

    seen = set()

    def add_seen(kmer):
        seen.add(kmer)
        #seen.add(revcomp(kmer))

    if exclude is not None:
        for pos in range(0, len(exclude) - K):
            add_seen(exclude[pos:pos + K - 1])

    seq = [random.choice('ACGT') for _ in range(K - 1)]  # do first K-1 bases
    add_seen(''.join(seq))

    while(len(seq) < length):
        next_base = random.choice('ACGT')
        next_kmer = ''.join(seq[-K + 2:] + [next_base])
        assert len(next_kmer) == K - 1
        if (next_kmer) not in seen:
            seq.append(next_base)
            add_seen(next_kmer)
        else:
            continue
    return ''.join(seq)


@pytest.fixture(params=list(range(500, 1600, 500)),
                ids=lambda val: '(L={0})'.format(val))
def random_sequence(request, K):

    def get(exclude=None):
        return get_random_sequence(request.param, K, exclude=exclude)

    return get


'''
# GRAPH STRUCTURE FIXTURE

These fixtures emit various graph structures with their corresponding
sequences and important nodes. They take a random sequence fixture and
a graph fixture, then consume sequence and generate k-mers accordingly.

We're using a bespoke but simple language to describe graph structures in the
docstrings of these tests. It is as follows:

    o: Node
    [x:y] <sequence>: Node at position in sequence, half-open intervals.
    (<name>), ([x:y] <name>): Named node, named node at position
    → : Edge
    ~~: Variable number of linear nodes.
'''

@pytest.fixture
def linear_subgraph(random_sequence):
    '''Sets up a simple linear path subgraph -- ie, a non-looping
    sequence with respect to K.

    sequence
    [0]→o→o~~o→o→[-1]
    '''
    sequence = random_sequence()

    return sequence


@pytest.fixture
def right_tip_subgraph(K, random_sequence):
    '''
    Sets up the sequences for a subgraph like so:

                             ([S-K+1:S+1]) tip
    sequence                ↗
    [0]→o→o~~o→([S-K:S] HDN)→o→o→o~~o→[-1]

    Where S is the mutation position.
    That is, it has a single branch at the S-Kth K-mer.

    Example:

    sequence: ACCTAATGGA
    tip:      ACCTAT
    S:             5
    '''
    sequence = random_sequence()
    S = random.randint(K+1, len(sequence)-K-1)
    tip_sequence = mutate_position(sequence[:S+1], -1)

    print(sequence, tip_sequence, sep='\n\n')
    assert tip_sequence[-2] == sequence[S-1]
    assert tip_sequence[-1] != sequence[S]

    return (sequence, tip_sequence), S


@pytest.fixture
def left_tip_subgraph(K, right_tip_subgraph):
    '''
    Sets up the sequences for a subgraph like so:

    tip ([0:K])
             ↘                    sequence
    [0]→o~~o→([S+1:S+K+1] HDN)→o→o~~o→[-1]

    Where S is the mutation position and S+1 is the starting
    position of the high degree (d=3) node.

    sequence: AGGTAATCCA
    tip:          TATCCA
    S:            4
    '''
    (sequence, tip), S = right_tip_subgraph
    sequence = sequence[::-1]
    tip = tip[::-1]
    S = len(sequence) - S - 1

    assert sequence[S+1] == tip[1]
    assert sequence[S] != tip[0]

    return (sequence, tip), S

@pytest.fixture
def snp_bubble_subgraph(request, right_tip_subgraph):
    '''
    Generate the sequences to create a SNP bubble.

                      ([S-K+1:S+1])→o~~o→([S:S+K]) wildtype
                    ↗                             ↘
    o~~([S-K:S] HDN)                               ([S+1:S+K+1] HDN)~~o
                    ↘                             ↗
                      ([S-K+1:S+1])→o~~o→([S:S+K]) mutant

    S is the location of the mutation relative to the shared start site of
    both sequences, which are labeled wildtype and mutant. 

    Example:

    wildtype: ATTCGTCAG
    mutant:   ATTGGTCAG
    S:           3
    '''

    (wildtype, mutant), S = right_tip_subgraph
    mutant += wildtype[S+1:]

    assert len(wildtype) == len(mutant)
    assert wildtype[S] != mutant[S]
    return (wildtype, mutant), S


