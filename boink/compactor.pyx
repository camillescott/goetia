# boink/compactor.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref

from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string

from boink.dbg cimport *
from boink.cdbg cimport *
from boink.prometheus cimport Instrumentation

from collections import namedtuple

Segment = namedtuple('Segment', ['sequence',
                                 'is_decision_kmer',
                                 'left_anchor',
                                 'right_anchor',
                                 'start',
                                 'length',
                                 'is_null'])


def display_segment_list(list segments):
    output = ''
    for i in range(len(segments)):
        cur_segment = segments[i]
        if cur_segment.is_null:
            output += '[x]'
        elif cur_segment.is_decision_kmer:
            output += '[D {0}]'.format(cur_segment.left_anchor)
        else:
            output += '[S {0} {1} {2}]'.format(cur_segment.left_anchor,
                                               cur_segment.right_anchor,
                                               cur_segment.length)
        if i != len(segments) - 1:
            output += '-'
    print(output)


include "compactor.tpl.pyx.pxi"

def make_streaming_compactor(dBG_Base graph, Instrumentation inst):
    return _make_streaming_compactor(graph, inst)
