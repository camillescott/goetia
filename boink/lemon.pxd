#
# Copyright 2011 Christoph Straehle, christoph.straehle@iwr.uni-heidelberg.de
# 
# Permission to use, modify and distribute this software is granted
#    * provided that this copyright notice appears in all copies. For
#    * precise terms see the accompanying LICENSE file.
#    *
#    * This software is provided "AS IS" with no warranty of any kind,
#    * express or implied, and with no claim as to its suitability for any
#    * purpose.
#
#############################################
# Part of this wrapper comes from the above author,
# and can be found here: https://github.com/cstraehl/cylemon/blob/master/cylemon/lemon/smart_graph.pxd
#############################################

from libcpp cimport bool
from libcpp.string cimport string


cdef extern from "<lemon/smart_graph.h>" namespace "lemon":
    cdef cppclass SmartDigraph


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass Arc
    cdef cppclass ArcIt
    cdef cppclass OutArcIt
    cdef cppclass InArcIt
    cdef cppclass Node
    cdef cppclass NodeIt


cdef extern from "<lemon/core.h>" namespace "lemon":
    cdef cppclass Invalid:
        Invalid()
    Invalid INVALID
    int countNodes[G](const G&)
    int countArcs[G](const G&)
    int countOutArcs[G](const G &, const Node &)
    int countInArcs[G](const G &, const Node &)



cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass Arc:
        inline Arc()
        inline Arc(Arc)
        inline Arc(ArcIt)
        inline Arc(OutArcIt)
        inline Arc(InArcIt)
        inline bint operator!=(Arc)
        inline bint operator!=(Invalid)
        inline bint operator<(Arc)
        inline bint operator==(Arc)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass ArcIt:
        inline Arc()
        inline Arc(Arc)
        inline ArcIt()
        inline ArcIt(ArcIt)
        inline ArcIt(SmartDigraph)
        inline ArcIt(SmartDigraph, Arc)
        inline ArcIt operator++()
        inline bint operator!=(Arc)
        inline bint operator!=(Invalid)
        inline bint operator<(Arc)
        inline bint operator==(Arc)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass OutArcIt:
        inline Arc()
        inline Arc(Arc)
        inline OutArcIt()
        inline OutArcIt(OutArcIt)
        inline OutArcIt(SmartDigraph&)
        inline OutArcIt(SmartDigraph&, Arc)
        inline OutArcIt(SmartDigraph&, Node)
        inline OutArcIt(SmartDigraph&, NodeIt)
        inline OutArcIt operator++()
        inline bint operator!=(Arc)
        inline bint operator!=(Invalid)
        inline bint operator<(Arc)
        inline bint operator==(Arc)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass InArcIt:
        inline Arc()
        inline Arc(Arc)
        inline InArcIt()
        inline InArcIt(InArcIt)
        inline InArcIt(SmartDigraph&)
        inline InArcIt(SmartDigraph&, Arc)
        inline InArcIt(SmartDigraph&, Node)
        inline InArcIt(SmartDigraph&, NodeItIt)
        inline InArcIt operator++()
        inline bint operator!=(Arc)
        inline bint operator!=(Invalid)
        inline bint operator<(Arc)
        inline bint operator==(Arc)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass Node:
        inline Node()
        inline Node(Node)
        inline Node(NodeIt)
        inline bint operator!=(Node)
        inline bint operator!=(Invalid)
        inline bint operator==(Invalid)
        inline bint operator!=(NodeIt)
        inline bint operator<(Node)
        inline bint operator==(Node)
        inline bint operator==(NodeIt)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass NodeIt:
        inline Node()
        inline Node(Node)
        inline NodeIt()
        inline NodeIt(SmartDigraph)
        inline NodeIt(SmartDigraph,Node)
        inline bint operator!=(Node)
        inline bint operator!=(Invalid)
        inline NodeIt operator++()
        inline bint operator<(Node)
        inline bint operator==(Node)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon":
    cdef cppclass SmartDigraph:                                         
        SmartDigraph()

        inline Node addNode()
        inline Arc addArc(Node, Node)

        void reserveNode(int)
        void reserveArc(int)
        void clear()

        # from Digraph concept:
        inline Node source(Arc)
        inline Node source(Arc&)
        inline Node source(ArcIt)
        inline Node source(OutArcIt)
        inline Node source(InArcIt)
        inline Node target(Arc)
        inline Node target(Arc&)
        inline Node target(ArcIt)
        inline Node target(OutArcIt)
        inline Node target(InArcIt)
        inline Node baseNode(OutArcIt)
        inline Node runningNode(OutArcIt)
        inline Node baseNode(InArcIt)
        inline Node runningNode(InArcIt)
        inline int id(Node)
        inline int id(NodeIt)
        inline int id(Arc)
        inline int id(ArcIt)
        inline int id(OutArcIt)
        inline int id(InArcIt)
        inline Node nodeFromId(int)
        inline Arc arcFromId(int)
        inline int maxNodeId()
        inline int maxArcId()
        inline Node oppositeNode(Node, Arc)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass NodeMap[T]:
        NodeMap(SmartDigraph)
        NodeMap(SmartDigraph, T)
        NodeMap(NodeMap[T])
        inline T& operator[](Node)
        inline T& operator[](NodeIt)
        inline void set(Node, T)
        inline void set(NodeIt, T)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass ArcMap[T]:
        ArcMap(SmartDigraph&)
        ArcMap(SmartDigraph&, T)
        inline T& operator[](Arc)
        inline T& operator[](ArcIt)
        inline T& operator[](OutArcIt)
        inline T& operator[](InArcIt)
        inline void set(Arc, T)
        inline void set(ArcIt, T)
        inline void set(InArcIt, T)
        inline void set(OutArcIt, T)


cdef extern from "<lemon/smart_graph.h>" namespace "lemon::SmartDigraph":
    cdef cppclass Snapshot:
        Snapshot()
        Snapshot(SmartDigraph&)
        void save(SmartDigraph&)
        void restore()


####################################
# End attributed code.
####################################


cdef extern from "lemon/maps.h" namespace "lemon":
    cdef cppclass CrossRefMap[G, K, V]:
        CrossRefMap(const G&)
        void set(const K&, const V&)
        const V operator[](const K&) const
        K get "operator()"(const V&) const

