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
    ctypedef pair[hash_t, hash_t] junction_t
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
        vector[junction_t] left_juncs
        vector[junction_t] right_juncs
        uint32_t count

        size_t degree() const
        size_t left_degree() const
        size_t right_degree() const

        bool has_left_junc(junction_t) const
        void add_left_junc(junction_t)
        bool has_right_junc(junction_t) const
        void add_right_junc(junction_t)

        string repr()

    cdef cppclass _UnitigNode "boink::UnitigNode" (_CompactNode):
        junction_t left_junc
        junction_t right_junc
        HashVector tags
        node_meta_t meta

        string repr()

        _UnitigNode(id_t, junction_t, junction_t, const string&)
        
    ctypedef _CompactNode * CompactNodePtr
    ctypedef _DecisionNode * DecisionNodePtr
    ctypedef _UnitigNode * UnitigNodePtr

    cdef cppclass _cDBG "boink::cDBG" (_KmerClient):
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

        _DecisionNode * get_dnode(hash_t)
        vector[_DecisionNode*] get_dnodes[ShifterType](const string&) except +ValueError

        _UnitigNode * get_unode(hash_t)
        _UnitigNode * get_unode(junction_t)
        _UnitigNode * get_unode_from_id(id_t)
        _DecisionNode * get_left_dnode(_UnitigNode*)
        _DecisionNode * get_right_dnode(_UnitigNode*)

        void write_adj_matrix(const string&) except +OSError
        void write_graphml(const string&) except +OSError

cdef extern from "boink/compactor.hh" namespace "boink" nogil:

    cdef cppclass _StreamingCompactor "boink::StreamingCompactor" [GraphType] (_AssemblerMixin[GraphType]):
        _cDBG cdbg

        _StreamingCompactor(GraphType *)

        #string compactify(const string&) except +ValueError
        #void compactify_right(Path&) 
        #void compactify_left(Path&)

        bool is_decision_kmer(uint8_t&)
        bool is_decision_kmer(const string&, uint8_t&) except +ValueError
        bool is_decision_kmer(const string&) except +ValueError
        void find_decision_kmers(const string&,
                                 vector[uint32_t]&,
                                 vector[hash_t]&,
                                 vector[NeighborBundle]&) except +ValueError
        void find_disturbed_dnodes(const string&,
                                   vector[_DecisionNode*]&,
                                   vector[NeighborBundle]&) except +ValueError

        bool insert_sequence(const string&,
                             vector[uint32_t]&,
                             vector[hash_t],
                             vector[NeighborBundle]&) except +ValueError

        void update_sequence(const string&) except +ValueError
        void update_cdbg(const string&) except +ValueError


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

