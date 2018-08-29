# cdbg.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

cimport cython
from libcpp.list cimport list as stdlist
from libcpp.memory cimport unique_ptr
from libcpp.pair cimport pair
from libcpp.unordered_set cimport unordered_set as uset
from libcpp.unordered_map cimport unordered_map as umap
from libcpp.set cimport set
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from boink.assembly cimport *
from boink.hashing cimport *

from boink.dbg cimport *
from boink.minimizers cimport _InteriorMinimizer
from boink.events cimport (_StreamingCompactorReport, _EventNotifier,
                           _EventListener, EventNotifier, EventListener)

cdef extern from "boink/cdbg.hh":
    cdef uint64_t NULL_ID
    cdef uint64_t UNITIG_START_ID

cdef extern from "boink/cdbg.hh" namespace "boink" nogil:

    ctypedef uint64_t id_t
    ctypedef pair[hash_t, hash_t] junction_t
    ctypedef vector[hash_t] HashVector

    ctypedef enum node_meta_t:
        FULL
        TIP
        ISLAND
        TRIVIAL

    ctypedef enum cDBGFormat:
        GRAPHML,
        EDGELIST,
        ADJMAT,
        FASTA,
        GFA1

    cdef const char * node_meta_repr(node_meta_t)

    cdef cppclass _CompactNode "boink::CompactNode":
        const id_t node_id
        string sequence

        _CompactNode(id_t, string)
        
        string revcomp()

        bool operator==(const _CompactNode&, const _CompactNode&)

    cdef cppclass _DecisionNode "boink::DecisionNode" (_CompactNode):
        uint32_t count()
        void incr_count()

        size_t degree() const
        size_t left_degree() const
        size_t right_degree() const

        string repr()

    cdef cppclass _UnitigNode "boink::UnitigNode" (_CompactNode):
        HashVector tags
        
        hash_t left_end()
        hash_t right_end()
        node_meta_t meta()

        string repr()

        _UnitigNode(id_t, junction_t, junction_t, const string&)
        
    ctypedef _CompactNode * CompactNodePtr
    ctypedef _DecisionNode * DecisionNodePtr
    ctypedef _UnitigNode * UnitigNodePtr

    cdef cppclass _cDBG "boink::cDBG" (_KmerClient, _EventListener):
        ctypedef umap[hash_t,unique_ptr[_DecisionNode]].const_iterator dnode_iter_t
        ctypedef umap[id_t,unique_ptr[_UnitigNode]].const_iterator unode_iter_t

        _cDBG(uint16_t K)

        uint64_t n_updates() const
        uint64_t n_unitig_nodes() const
        uint64_t n_decision_nodes() const
        uint64_t n_tags() const

        unode_iter_t unodes_begin() const
        unode_iter_t unodes_end() const
        dnode_iter_t dnodes_begin() const
        dnode_iter_t dnodes_end() const

        _DecisionNode * query_dnode(hash_t)
        bool query_dnode_marker(hash_t)
        vector[_DecisionNode*] query_dnodes[ShifterType](const string&) except +ValueError

        _UnitigNode * query_unode_tag(hash_t)
        _UnitigNode * query_unode_end(hash_t)
        _UnitigNode * query_unode_id(id_t)
    
        void write(const string&, cDBGFormat) except +OSError
        void write_adj_matrix(const string&) except +OSError
        void write_graphml(const string&) except +OSError

cdef extern from "boink/compactor.hh" namespace "boink" nogil:

    cdef struct compact_segment:
        hash_t left_anchor
        hash_t right_anchor
        size_t start_pos
        size_t length
        bool is_decision_kmer

        compact_segment()
        const bool is_null() 

    cdef cppclass _StreamingCompactor "boink::StreamingCompactor" [GraphType] (_AssemblerMixin[GraphType], _EventNotifier):
        _cDBG cdbg

        _StreamingCompactor(GraphType *)

        #string compactify(const string&) except +ValueError
        #void compactify_right(Path&) 
        #void compactify_left(Path&)

        void wait_on_updates()

        bool is_decision_kmer(uint8_t&)
        bool is_decision_kmer(const string&, uint8_t&) except +ValueError
        bool is_decision_kmer(const string&) except +ValueError

        void find_decision_kmers(const string&,
                                 vector[uint32_t]&,
                                 vector[hash_t]&,
                                 vector[NeighborBundle]&) except +ValueError

        void update_sequence(const string&) except +ValueError

        void find_new_segments(const string&, # sequence to add
                               set[hash_t]&, # all new k-mers
                               deque[compact_segment]&, # new segments
                               set[hash_t]&, # new decision k-mers
                               deque[NeighborBundle]& # decision neighbors
                               ) except +ValueError

        _StreamingCompactorReport* get_report()

cdef cDBGFormat convert_format(str graph_format) except *


cdef class CompactNode:
    cdef _CompactNode * _cn_this

    @staticmethod
    cdef CompactNode _wrap(_CompactNode *)


cdef class DecisionNode(CompactNode):
    cdef _DecisionNode * _dn_this

    @staticmethod
    cdef DecisionNode _wrap(_DecisionNode *)


cdef class UnitigNode(CompactNode):
    cdef _UnitigNode * _un_this

    @staticmethod
    cdef UnitigNode _wrap(_UnitigNode *)


cdef class cDBG:
    cdef _cDBG * _this

    @staticmethod
    cdef cDBG _wrap(_cDBG *)


cdef class StreamingCompactor:
    cdef DefaultDBG * _graph
    # TODO jinja template StreamingCompactor template args
    cdef unique_ptr[_StreamingCompactor[DefaultDBG]] _sc_this
    cdef public cDBG cdbg
    cdef public EventNotifier Notifier

