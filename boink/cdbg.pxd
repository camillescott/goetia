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
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from boink.assembly cimport *
from boink.hashing cimport *

from boink.dbg cimport *
from boink.minimizers cimport _InteriorMinimizer

cdef extern from "boink/cdbg.hh":
    cdef uint64_t NULL_ID
    cdef uint64_t UNITIG_START_ID

cdef extern from "boink/cdbg.hh" namespace "boink" nogil:

    ctypedef uint64_t id_t
    ctypedef vector[hash_t] HashVector

    ctypedef enum node_meta_t:
        FULL
        TIP
        ISLAND
        TRIVIAL

    cdef const char * node_meta_repr(node_meta_t)

    cdef cppclass _CompactNode "boink::CompactNode":
        const id_t node_id
        string sequence

        _CompactNode(id_t, string)
        
        string revcomp()

        bool operator==(const _CompactNode&, const _CompactNode&)

    cdef cppclass _DecisionNode "boink::DecisionNode" (_CompactNode):
        vector[id_t] in_edges
        vector[id_t] out_edges
        uint32_t count

        size_t degree() const
        size_t in_degree() const
        size_t out_degree() const

        bool has_in_edge(id_t) const
        void remove_in_edge(id_t)
        bool has_out_edge(id_t) const
        void remove_out_edge(id_t)

        string repr()

    cdef cppclass _UnitigNode "boink::UnitigNode" (_CompactNode):
        _DecisionNode * left_dnode
        _DecisionNode * right_dnode
        hash_t left_hash
        hash_t right_hash
        HashVector tags

        string repr()

        _UnitigNode(id_t, hash_t, hash_t, const string&)
        
    ctypedef _CompactNode * CompactNodePtr
    ctypedef _DecisionNode * DecisionNodePtr
    ctypedef _UnitigNode * UnitigNodePtr

    cdef cppclass _cDBG "boink::cDBG" (_KmerClient):
        _cDBG(uint16_t K)

        uint64_t n_updates() const
        uint64_t n_unitig_nodes() const
        uint64_t n_decision_nodes() const

        _DecisionNode * get_decision_node(hash_t)
        vector[_DecisionNode*] get_decision_nodes[ShifterType](const string&) except +ValueError

        _UnitigNode * get_unitig_node(hash_t)
        _UnitigNode * get_unitig_node_from_id(id_t)

    cdef cppclass _StreamingCompactor "boink::StreamingCompactor" [GraphType] (_AssemblerMixin[GraphType]):
        _cDBG cdbg

        _StreamingCompactor(GraphType *)

        #string compactify(const string&) except +ValueError
        #void compactify_right(Path&) 
        #void compactify_left(Path&)

        bool is_decision_node(uint8_t&)
        bool is_decision_node(const string&, uint8_t&) except +ValueError
        bool is_decision_node(const string&) except +ValueError
        void find_decision_nodes(const string&,
                                 vector[uint32_t]&,
                                 vector[hash_t]&,
                                 vector[NeighborBundle]&) except +ValueError
        void find_disturbed_decision_nodes(const string&,
                                           vector[_DecisionNode*]&,
                                           vector[NeighborBundle]&) except +ValueError

        bool insert_sequence(const string&,
                             vector[uint32_t]&,
                             vector[hash_t],
                             vector[NeighborBundle]&) except +ValueError

        void update(const string&) except +ValueError


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

