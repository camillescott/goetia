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
from libcpp.set cimport set
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from boink.assembly cimport *
from boink.hashing cimport *

from boink.dbg cimport *
from boink.minimizers cimport _InteriorMinimizer
from boink.events cimport (_StreamingCompactorReport, _EventNotifier,
                           _EventListener, EventNotifier, EventListener)

from boink.sparsepp cimport sparse_hash_set, sparse_hash_map

cdef extern from "boink/cdbg.hh":
    cdef uint64_t NULL_ID
    cdef uint64_t UNITIG_START_ID

cdef extern from "boink/cdbg.hh" namespace "boink" nogil:

    ctypedef uint64_t id_t
    ctypedef pair[hash_t, hash_t] junction_t
    ctypedef vector[hash_t] HashVector

    ctypedef enum node_meta_t:
        FULL,
        TIP,
        ISLAND,
        CIRCULAR,
        LOOP,
        TRIVIAL,
        DECISION

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

        _CompactNode(id_t, string, node_meta_t)
        
        string revcomp()
        node_meta_t meta()

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

        string repr()

        _UnitigNode(id_t, junction_t, junction_t, const string&)
        
    ctypedef _CompactNode * CompactNodePtr
    ctypedef _DecisionNode * DecisionNodePtr
    ctypedef _UnitigNode * UnitigNodePtr
    
    cdef cppclass _cDBG "boink::cDBG" [GraphType] (_KmerClient, _EventNotifier):
        ctypedef sparse_hash_map[hash_t,unique_ptr[_DecisionNode]].const_iterator dnode_iter_t
        ctypedef sparse_hash_map[id_t,unique_ptr[_UnitigNode]].const_iterator unode_iter_t

        _cDBG(uint16_t K)

        GraphType * get_dbg()

        uint64_t n_updates() const
        uint64_t n_unitig_nodes() const
        uint64_t n_decision_nodes() const
        uint64_t n_tags() const
        uint64_t n_unitig_ends() const

        unode_iter_t unodes_begin() const
        unode_iter_t unodes_end() const
        dnode_iter_t dnodes_begin() const
        dnode_iter_t dnodes_end() const

        _DecisionNode * query_dnode(hash_t)
        bool has_dnode(hash_t)
        vector[_DecisionNode*] query_dnodes(const string&) except +ValueError
        pair[vector[CompactNodePtr], vector[CompactNodePtr]] find_dnode_neighbors(_DecisionNode *)

        CompactNodePtr query_cnode(hash_t)

        _UnitigNode * query_unode_tag(hash_t)
        _UnitigNode * query_unode_end(hash_t)
        _UnitigNode * query_unode_id(id_t)
        pair[DecisionNodePtr, DecisionNodePtr] find_unode_neighbors(_UnitigNode*)
    
        void validate(const string&) except+ OSError
        void write(const string&, cDBGFormat) except +OSError
        void write_adj_matrix(const string&) except +OSError
        void write_graphml(const string&) except +OSError

    cdef cppclass _AsyncCDBG "boink::AsyncCDBG" [GraphType] (_cDBG, _EventListener):
        _AsyncCDBG(uint16_t K)


cdef cDBGFormat convert_format(str graph_format) except *


cdef class CompactNodeView:
    cdef _CompactNode * _cn_this

    @staticmethod
    cdef CompactNodeView _wrap(_CompactNode *)


cdef class DecisionNodeView(CompactNodeView):
    cdef _DecisionNode * _dn_this

    @staticmethod
    cdef DecisionNodeView _wrap(_DecisionNode *)


cdef class UnitigNodeView(CompactNodeView):
    cdef _UnitigNode * _un_this

    @staticmethod
    cdef UnitigNodeView _wrap(_UnitigNode *)


cdef class CompactNode:
    cdef id_t node_id
    cdef string sequence


cdef class DecisionNode(CompactNode):
    cdef uint32_t count

    @staticmethod
    cdef DecisionNode _create(DecisionNodeView view)


cdef class UnitigNode(CompactNode):
    cdef hash_t left_end
    cdef hash_t right_end
    cdef node_meta_t meta

    @staticmethod
    cdef UnitigNode _create(UnitigNodeView view)

include "cdbg.tpl.pxd.pxi"
