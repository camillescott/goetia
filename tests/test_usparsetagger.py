# boink/tests/test_usparsetagger.py
# Copyright (C) 2020 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import itertools
import pytest

from tests.utils import *

from boink import libboink
from boink.dbg import dBG
from boink.hashing import CanUnikmerShifter
from boink.storage import SparseppSetStorage


@pytest.fixture
def tagger(ksize, store):
    tagger_type = libboink.cdbg.USparseGraph[type(store)].Graph
    hasher = CanUnikmerShifter.build(ksize, 7)
    graph = dBG[type(store), CanUnikmerShifter].build(store, hasher)
    tagger = tagger_type.build(graph, hasher.ukhs_map)
    return tagger


@using(ksize=31, length=100)
@exact_backends()
class TestFindNewExtensions:

    def test_new_whole_sequence(self, ksize, length, tagger, linear_path):
        sequence = linear_path()
        extensions = tagger.find_new_extensions(sequence)

        assert len(extensions) == len(sequence) - ksize + 1
        for h, pos, (left, right) in extensions:
            assert len(left) == 4
            assert len(right) == 4
    
    def test_new_except_front(self, ksize, length, tagger, linear_path):
        sequence = linear_path()
        tagger.dbg.insert(sequence[:ksize])
        extensions = tagger.find_new_extensions(sequence)

        assert len(extensions) == len(sequence) - ksize

    def test_new_single_kmer(self, ksize, length, tagger, linear_path):
        sequence = linear_path()
        position = 10

        tagger.dbg.insert_sequence(sequence[:position + ksize - 1])
        tagger.dbg.insert_sequence(sequence[position+1:])

        extensions = tagger.find_new_extensions(sequence)

        assert len(extensions) == 1
        h, pos, (left, right) = extensions[0]
        assert h == tagger.dbg.hash(sequence[position:position+ksize])
        assert pos == position


@using(ksize=31, length=100)
@exact_backends()
class TestFilterNewExtensions:
    
    def test_new_whole_sequence(self, ksize, length, tagger, linear_path):
        ''' Test that a sequence consisting of new k-mers produces a properly filtered
        chain of neighbors.
        '''
        sequence = linear_path()
        extensions = tagger.filter_new_extensions(tagger.find_new_extensions(sequence))

        assert len(extensions) == len(sequence) - ksize + 1
        for h, pos, (left, right) in extensions[1:-1]:
            assert len(left) == 1
            assert len(right) == 1
        
        fh, fpos, (fleft, fright) = extensions[0]
        assert len(fleft) == 0
        assert len(fright) == 1

        bh, bpos, (bleft, bright) = extensions[-1]
        assert len(bleft) == 1
        assert len(bright) == 0

    def test_new_except_front(self, ksize, length, tagger, linear_path):
        ''' Test that a sequence where a k-mer is non-new returns a properly
        filtered chain of neighbors.
        '''
        sequence = linear_path()
        tagger.dbg.insert(sequence[:ksize])
        extensions = tagger.filter_new_extensions(tagger.find_new_extensions(sequence))

        assert len(extensions) == len(sequence) - ksize

        # the first k-mer exists in the graph, so all the new k-mers
        # except for the last should have degree (1,1)
        for h, pos, (left, right) in extensions[:-1]:
            assert len(left) == 1
            assert len(right) == 1

        bh, bpos, (bleft, bright) = extensions[-1]
        assert len(bleft) == 1
        assert len(bright) == 0


@using(ksize=31, length=100)
@exact_backends()
class TestBuildNewSegments:

    def test_new_whole_sequence(self, ksize, length, tagger, linear_path):
        ''' Test that a sequence consisting entirely of new k-mers produces
        one segment chain matching that produced by TestFilteredNewExtensions::test_new_whole_sequence.
        '''
        sequence = linear_path()
        segments = tagger.build_new_segments(
                       tagger.filter_new_extensions(
                           tagger.find_new_extensions(sequence)
                       )
                   )

        assert len(segments) == 1

        extensions = segments[0]
        assert len(extensions) == len(sequence) - ksize + 1
        for h, pos, (left, right) in extensions[1:-1]:
            assert len(left) == 1
            assert len(right) == 1
        
        fh, fpos, (fleft, fright) = extensions[0]
        assert len(fleft) == 0
        assert len(fright) == 1

        bh, bpos, (bleft, bright) = extensions[-1]
        assert len(bleft) == 1
        assert len(bright) == 0

    def test_new_split_sequence(self, ksize, length, tagger, linear_path):
        ''' Test that a sequence with a non-new k-mer in the center returns
        two segments.
        '''
        sequence = linear_path()
        position = 10

        tagger.dbg.insert(sequence[position:position+ksize])
        segments = tagger.build_new_segments(
                       tagger.filter_new_extensions(
                           tagger.find_new_extensions(sequence)
                       )
                   )

        assert len(segments) == 2
        assert len(segments[0]) == 10
        assert len(segments[1]) == len(sequence) - ksize - 10

        lsegment, rsegment = segments

        _, _, (_, lsegment_right) = lsegment[-1]
        assert len(lsegment_right) == 1
        assert lsegment_right[0].hash == tagger.dbg.hash(sequence[position:position+ksize])

        _, _, (rsegment_left, _) = rsegment[0]
        assert len(rsegment_left) == 1
        assert rsegment_left[0].hash == tagger.dbg.hash(sequence[position:position+ksize])

    def test_new_single_kmer(self, ksize, length, tagger, linear_path):
        sequence = linear_path()
        position = 10

        tagger.dbg.insert_sequence(sequence[:position + ksize - 1])
        tagger.dbg.insert_sequence(sequence[position+1:])

        segments = tagger.build_new_segments(
                       tagger.filter_new_extensions(
                           tagger.find_new_extensions(sequence)
                       )
                   )

        assert len(segments) == 1
        assert len(segments[0]) == 1
        
        h, pos, (left, right) = segments[0][0]
        assert h == tagger.dbg.hash(sequence[position:position+ksize])
        assert pos == position
        assert len(left) == 1
        assert len(right) == 1


@using(ksize=31, length=100)
@exact_backends()
class TestSplitNewSegments:

    def test_new_whole_sequence(self, ksize, length, tagger, linear_path):
        ''' Test that a sequence consisting entirely of new k-mers produces
        one segment chain matching that produced by TestBuildNewSegments::test_new_whole_sequence.
        '''
        sequence = linear_path()
        segments = tagger.split_new_segments(
                       tagger.build_new_segments(
                           tagger.filter_new_extensions(
                               tagger.find_new_extensions(sequence)
                           )
                       )
                   )   

        assert len(segments) == 1

        extensions = segments[0]
        assert len(extensions) == len(sequence) - ksize + 1
        for h, pos, (left, right) in extensions[1:-1]:
            assert len(left) == 1
            assert len(right) == 1
        
        fh, fpos, (fleft, fright) = extensions[0]
        assert len(fleft) == 0
        assert len(fright) == 1

        bh, bpos, (bleft, bright) = extensions[-1]
        assert len(bleft) == 1
        assert len(bright) == 0

    def test_fwd_decision_split(self, ksize, length, tagger, right_fork, check_fp):
        ''' Test that segments get split on FWD decision k-mers.
        '''
        (core, branch), pivot = right_fork()
        check_fp()

        # insert the branch: it does *not* contain the decision k-mer
        tagger.dbg.insert_sequence(branch)
        # insert core, which does contain the decision k-mer
        segments = tagger.split_new_segments(
                       tagger.build_new_segments(
                           tagger.filter_new_extensions(
                               tagger.find_new_extensions(core)
                           )
                       )
                   )   

        print([(pos, len(left), len(right)) for h, pos, (left, right) in segments[0]])
        assert len(segments) == 2

        lsegment, rsegment = segments
        lsegment_rflank, _, (lsegment_rflank_left, lsegment_rflank_right) = lsegment[-1]
        assert len(lsegment_rflank_right) == 2
        assert lsegment_rflank == tagger.dbg.hash(core[pivot:pivot+ksize])

        assert rsegment[0][0] == tagger.dbg.hash(core[pivot+1:pivot+ksize+1])

    def test_rev_decision_split(self, ksize, length, tagger, left_fork, check_fp):
        ''' Test that segments get split on FWD decision k-mers.
        '''
        (core, branch), pivot = left_fork()
        check_fp()

        # insert the branch: it does *not* contain the decision k-mer
        tagger.dbg.insert_sequence(branch)
        # insert core, which does contain the decision k-mer
        segments = tagger.split_new_segments(
                       tagger.build_new_segments(
                           tagger.filter_new_extensions(
                               tagger.find_new_extensions(core)
                           )
                       )
                   )   

        print([(pos, len(left), len(right)) for h, pos, (left, right) in segments[0]])
        assert len(segments) == 2

        lsegment, rsegment = segments
        rsegment_lflank, _, (rsegment_lflank_left, rsegment_lflank_right) = rsegment[0]
        lsegment_rflank, _, (lsegment_rflank_left, lsegment_rflank_right) = lsegment[-1]

        assert len(lsegment_rflank_right) == 1
        assert lsegment_rflank == tagger.dbg.hash(core[pivot-1:pivot+ksize-1])

        assert len(rsegment_lflank_left) == 2
        assert rsegment_lflank == tagger.dbg.hash(core[pivot:pivot+ksize])

    def test_full_decision_split(self, ksize, length, tagger, full_decision, check_fp):
        ''' Test the a full decision k-mer (in and out-degree > 2) gets its own segment.
        '''

        core, (lbranches, rbranches) = full_decision()
        #check_fp()

        for seq in itertools.chain(lbranches, rbranches):
            tagger.dbg.insert_sequence(seq)
        
        segments = tagger.split_new_segments(
                       tagger.build_new_segments(
                           tagger.filter_new_extensions(
                               tagger.find_new_extensions(core)
                           )
                       )
                   ) 

        print(len(lbranches))
        print([(pos, len(left), len(right)) for h, pos, (left, right) in segments[0]])
        print([(pos, len(left), len(right)) for h, pos, (left, right) in segments[1]])
        assert len(segments) == 3