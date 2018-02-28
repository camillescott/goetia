import pytest

from boink.dbg import make_dBG
from boink.tests.utils import *
from khmer.tests.graph_structure_fixtures import *


@pytest.fixture(params=['BitStorage', 'ByteStorage', 'NibbleStorage'])
def dbg_type(request, ksize):

    def build(*args):
        if not args:
            starting_size, n_tables = 1000000, 4
        else:
            starting_size, n_tables = args
        return make_dBG(ksize, starting_size, n_tables, storage=request.param)

    return build


@using_ksize([21, 51, 101])
def test_presence(dbg_type, ksize, random_sequence):
    # basic get/add test
    tt = dbg_type()

    for kmer in kmers(random_sequence(), ksize):

        hashval = tt.hash(kmer)

        assert tt.get(kmer) == 0
        assert tt.get(hashval) == 0

        tt.add(kmer)
        assert tt.get(kmer) == 1
        assert tt.get(hashval) == 1

        tt.add(kmer)
        # Node* types can only tell presence/absence
        if 'Bit' in tt.__class__.__name__:
            assert tt.get(kmer) == 1
            assert tt.get(hashval) == 1
        else:
            assert tt.get(kmer) == 2
            assert tt.get(hashval) == 2


@using_ksize([21,151])
def test_n_occupied(dbg_type, ksize):
    # basic get/add test
    tt = dbg_type()

    kmer = 'G' * ksize

    assert tt.n_occupied == 0
    assert tt.n_unique == 0

    tt.add(kmer)
    assert tt.n_occupied == 1
    assert tt.n_unique == 1

    tt.add(kmer)
    # the CQF implementation we use can use more than one slot to represent
    # counts for a single kmer
    if not "QF" in tt.__class__.__name__:
        assert tt.n_occupied == 1
    else:
        assert tt.n_occupied == 2
    assert tt.n_unique == 1


'''
def test_bad_create(Tabletype):
    # creation should fail w/bad parameters
    try:
        tt = Tabletype(5, [])
    except ValueError as err:
        assert 'tablesizes needs to be one or more numbers' in str(err)
'''

@using_ksize([21,51,81])
def test_get_ksize(dbg_type, ksize):
    # ksize() function.
    kh = dbg_type()
    assert kh.K == ksize


def test_hash(dbg_type, ksize):
    # hashing of strings -> numbers.
    kh = dbg_type()
    x = kh.hash("ATGGC")
    assert type(x) == int


def test_hash_bad_dna(dbg_type, ksize):
    # hashing of bad dna -> succeeds w/o complaint
    kh = dbg_type()

    x = kh.hash("ATGYC")


def test_hash_bad_length(dbg_type, ksize):
    # hashing of too long should ignore extra sequence
    kh = dbg_type()
    test_kmer = 'A' * ksize
    assert kh.hash(test_kmer) == kh.hash(test_kmer + 'TTTT')


# TODO add test for k-mer too short

'''

def test_hashsizes(AnyTabletype):
    # hashsizes method.
    kh = AnyTabletype(5)
    assert (kh.hashsizes() == PRIMES_1m or
            # CQF allocates some extra slots beyond what you request
            # exactly how many extra is an implementation detail
            kh.hashsizes()[0] >= QF_SIZE)
'''


@using_ksize(5)
def test_add_hashval(dbg_type, ksize):
    # test add(hashval)
    kh = dbg_type()
    x = kh.hash("ATGGC")
    y = kh.add(x)
    assert y

    z = kh.get(x)
    assert z == 1


@using_ksize(5)
def test_add_dna_kmer(dbg_type, ksize):
    # test add(dna)
    kh = dbg_type()
    x = kh.add("ATGGC")
    assert x

    z = kh.get("ATGGC")
    assert z == 1


@using_ksize(5)
def test_get_hashval(dbg_type, ksize):
    # test get(hashval)
    kh = dbg_type()
    hashval = kh.hash("ATGGC")
    kh.add(hashval)

    z = kh.get(hashval)
    assert z == 1


@using_ksize(5)
def test_get_hashval_rc(dbg_type, ksize):
    # fw and rc should NOT be the same on this table
    kh = dbg_type()
    hashval = kh.hash("ATGC")
    rc = kh.hash("GCAT")

    assert hashval != rc

'''
def test_get_dna_kmer(AnyTabletype):
    # test get(dna)
    kh = AnyTabletype(5)
    hashval = kh.hash("ATGGC")
    kh.add(hashval)

    z = kh.get("ATGGC")
    assert z == 1


def test_get_bad_dna_kmer(AnyTabletype):
    # test get(dna) with bad dna; should be fine.
    kh = AnyTabletype(5)

    kh.hash("ATYGC")


def test_consume_and_count(AnyTabletype):
    tt = AnyTabletype(6)

    x = "ATGCCGATGCA"
    num_kmers = tt.consume(x)
    assert num_kmers == len(x) - tt.ksize() + 1   # num k-mers consumed

    for start in range(len(x) - 6 + 1):
        assert tt.get(x[start:start + 6]) == 1


def test_consume_and_count_bad_dna(AnyTabletype):
    # while we don't specifically handle bad DNA, we should at least be
    # consistent...
    tt = AnyTabletype(6)

    x = "ATGCCGNTGCA"
    num_kmers = tt.consume(x)

    for start in range(len(x) - 6 + 1):
        assert tt.get(x[start:start + 6]) == 1


def test_consume_short(AnyTabletype):
    # raise error on too short when consume is run
    tt = AnyTabletype(6)

    x = "ATGCA"
    with pytest.raises(ValueError):
        tt.consume(x)


def test_get_kmer_counts(AnyTabletype):
    hi = AnyTabletype(6)

    hi.consume("AAAAAA")
    counts = hi.get_kmer_counts("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] == 1

    hi.consume("AAAAAA")
    counts = hi.get_kmer_counts("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] >= 1

    hi.consume("AAAAAT")
    counts = hi.get_kmer_counts("AAAAAAT")
    print(counts)
    assert len(counts) == 2
    assert counts[0] >= 1
    assert counts[1] == 1


def test_get_kmer_hashes(AnyTabletype):
    hi = AnyTabletype(6)

    hashes = hi.get_kmer_hashes("ACGTGCGT")
    print(hashes)
    assert len(hashes) == 3
    assert hashes[0] == hi.hash("ACGTGC")
    assert hashes[1] == hi.hash("CGTGCG")
    assert hashes[2] == hi.hash("GTGCGT")


def test_get_min_count(AnyTabletype):
    hi = AnyTabletype(6)

    # master string, 3 k-mers
    x = "ACGTGCGT"

    hi.add("ACGTGC")  # 3
    hi.add("ACGTGC")
    hi.add("ACGTGC")

    hi.add("CGTGCG")  # 1

    hi.add("GTGCGT")  # 2
    hi.add("GTGCGT")

    counts = hi.get_kmer_counts(x)
    assert hi.get_min_count(x) == min(counts)
    assert hi.get_max_count(x) == max(counts)
    med, _, _ = hi.get_median_count(x)
    assert med == list(sorted(counts))[len(counts) // 2]


def test_get_kmers(AnyTabletype):
    hi = AnyTabletype(6)

    kmers = hi.get_kmers("AAAAAA")
    assert kmers == ["AAAAAA"]

    kmers = hi.get_kmers("AAAAAAT")
    assert kmers == ["AAAAAA", "AAAAAT"]

    kmers = hi.get_kmers("AGCTTTTC")
    assert kmers == ['AGCTTT', 'GCTTTT', 'CTTTTC']


def test_trim_on_abundance(AnyTabletype):
    hi = AnyTabletype(6)

    x = "ATGGCAGTAGCAGTGAGC"
    hi.consume(x[:10])

    (y, pos) = hi.trim_on_abundance(x, 1)
    assert pos == 10
    assert x[:pos] == y


def test_trim_below_abundance(AnyTabletype):
    hi = AnyTabletype(6)

    x = "ATGGCAGTAGCAGTGAGC"
    x_rc = screed.rc(x)
    hi.consume(x_rc[:10])

    print(len(x))

    (y, pos) = hi.trim_below_abundance(x, 0)
    assert pos == len(x) - hi.ksize() + 1
    assert x[:pos] == y


DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"


def test_find_spectral_error_positions(AnyTabletype):
    kh = AnyTabletype(8)

    kh.consume(DNA[:30])

    for n in range(len(DNA) - 8 + 1):
        print(n, kh.get(DNA[n:n + 8]))

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [30], posns


def test_find_spectral_error_positions_6(AnyTabletype):
    kh = AnyTabletype(8)

    kh.consume(DNA[1:])

    for n in range(len(DNA) - 8 + 1):
        print(n, kh.get(DNA[n:n + 8]))

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [0], posns


def test_find_spectral_error_positions_5(AnyTabletype):
    kh = AnyTabletype(8)

    kh.consume(DNA[:10])
    kh.consume(DNA[11:])

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [10], posns


def test_consume_seqfile_reads_parser(AnyTabletype):
    kh = AnyTabletype(5)
    rparser = ReadParser(utils.get_test_data('test-fastq-reads.fq'))

    kh.consume_seqfile(rparser)

    kh2 = AnyTabletype(5)
    for record in screed.open(utils.get_test_data('test-fastq-reads.fq')):
        kh2.consume(record.sequence)

    assert kh.get('CCGGC') == kh2.get('CCGGC')


def test_consume_seqfile(AnyTabletype):
    kh = AnyTabletype(5)
    kh.consume_seqfile(utils.get_test_data('test-fastq-reads.fq'))

    kh2 = AnyTabletype(5)
    for record in screed.open(utils.get_test_data('test-fastq-reads.fq')):
        kh2.consume(record.sequence)

    assert kh.get('CCGGC') == kh2.get('CCGGC')


def test_save_load(Tabletype):
    kh = Tabletype(5)
    ttype = type(kh)
    savefile = utils.get_temp_filename('tablesave.out')

    # test add(dna)
    x = kh.add("ATGGC")
    z = kh.get("ATGGC")
    assert z == 1

    kh.save(savefile)

    # should we provide a single load function here? yes, probably. @CTB
    loaded = ttype.load(savefile)

    z = loaded.get('ATGGC')
    assert z == 1


def test_get_bigcount(Tabletype):
    # get_bigcount should return false by default
    tt = Tabletype(12)

    assert not tt.get_use_bigcount()


def test_set_bigcount(Tabletype):
    supports_bigcount = [Countgraph, Counttable, CyclicCounttable]
    tt = Tabletype(12)

    if type(tt) in supports_bigcount:
        tt.set_use_bigcount(True)

        for i in range(300):
            tt.add('G' * 12)
        assert tt.get('G' * 12) == 300

    else:
        with pytest.raises(ValueError):
            tt.set_use_bigcount(True)


def test_abund_dist_A(AnyTabletype):
    A_filename = utils.get_test_data('all-A.fa')

    kh = AnyTabletype(4)
    tracking = Nodegraph(4, 1, 1, primes=PRIMES_1m)

    kh.consume_seqfile(A_filename)
    dist = kh.abundance_distribution(A_filename, tracking)

    print(dist[:10])
    assert sum(dist) == 1
    assert dist[0] == 0


def test_abund_dist_A_readparser(AnyTabletype):
    A_filename = utils.get_test_data('all-A.fa')
    rparser = ReadParser(A_filename)

    kh = AnyTabletype(4)
    tracking = Nodegraph(4, 1, 1, primes=PRIMES_1m)

    kh.consume_seqfile(A_filename)
    dist = kh.abundance_distribution(rparser, tracking)

    print(dist[:10])
    assert sum(dist) == 1
    assert dist[0] == 0
'''
