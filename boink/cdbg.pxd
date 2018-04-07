# cdbg.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

cimport cython
from libcpp.memory cimport shared_ptr
from libcpp.list cimport list as stdlist
from libcpp.pair cimport pair
from libcpp.unordered_set cimport unordered_set as uset
from libcpp.unordered_map cimport unordered_map as umap
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from boink.assembly cimport *
from boink.hashing cimport *
from boink.dbg cimport *

cdef extern from "boink/cdbg.hh":
    cdef uint64_t NULL_ID
    cdef uint64_t UNITIG_START_ID

cdef extern from "boink/cdbg.hh" namespace "boink" nogil:

    ctypedef uint64_t id_t
    ctypedef pair[hash_t, id_t] HashIDPair
    ctypedef uset[hash_t] UHashSet
    ctypedef vector[hash_t] HashVector
    ctypedef umap[hash_t, id_t] HashIDMap
    ctypedef uset[id_t] IDSet

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

    cdef cppclass _UnitigNode "boink::UnitigNode" (_CompactNode):
        id_t in_id
        id_t out_id
        HashVector tags

        _UnitigNode(id_t, string, uint64_t)

    cdef cppclass _DecisionNode "boink::DecisionNode" (_CompactNode):
        vector[id_t] in_edges
        vector[id_t] out_edges
        uint32_t count

        uint8_t degree() const
        uint8_t in_degree() const
        uint8_t out_degree() const
        

    ctypedef vector[_CompactNode] CompactNodeVector
    ctypedef shared_ptr[_CompactNode] CompactNodePtr
    ctypedef shared_ptr[_DecisionNode] DecisionNodePtr
    ctypedef shared_ptr[_UnitigNode] UnitigNodePtr

    cdef cppclass _cDBG "boink::cDBG" (_KmerClient):
        _cDBG(uint16_t K)

    cdef cppclass _StreamingCompactor "boink::StreamingCompactor" [GraphType] (_AssemblerMixin[GraphType]):

        _StreamingCompactor(shared_ptr[GraphType])

        string compactify(const string&) except +ValueError
        void compactify_right(Path&) 
        void compactify_left(Path&)

        bool is_decision_node(uint8_t&)
        bool is_decision_node(const string&, uint8_t&) except +ValueError
        bool is_decision_node(const string&) except +ValueError

        bool insert_sequence(const string&,
                             vector[uint32_t]&,
                             vector[hash_t],
                             vector[NeighborBundle]&) except +ValueError

        void update(const string&) except +ValueError



cdef class CompactNode:
    cdef CompactNodePtr _cn_this

    @staticmethod
    cdef CompactNode _wrap(CompactNodePtr)


cdef class DecisionNode(CompactNode):
    cdef DecisionNodePtr _dn_this

    @staticmethod
    cdef DecisionNode _wrap(DecisionNodePtr)


cdef class UnitigNode(CompactNode):
    cdef UnitigNodePtr _un_this

    @staticmethod
    cdef UnitigNode _wrap(UnitigNodePtr)


cdef class StreamingCompactor:

    cdef shared_ptr[DefaultDBG] _graph
    # TODO jinja template StreamingCompactor template args
    cdef shared_ptr[_StreamingCompactor[DefaultDBG]] _sc_this

