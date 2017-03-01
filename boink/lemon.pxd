from libcpp cimport bool
from libcpp.string cimport string


cdef extern from "lemon/core.h" namespace "lemon":

    cdef cppclass Invalid "lemon::Invalid":
        bool operator==(Invalid)
        bool operator!=(Invalid)
        bool operator<(Invalid)
    cdef Invalid INVALID "lemon::INVALID"
        

cdef extern from "lemon/smart_graph.h" namespace "lemon":
    cdef cppclass SmartDigraph:
        SmartDigraph()

        cppclass Node:
            Node()
            Node(const Node &)
            Node(Invalid)
            bool operator==(Node) const
            bool operator!=(Node) const
            bool operator<(Node) const

        cppclass NodeIt(Node):
            NodeIt()
            NodeIt(Invalid)
            NodeIt(const SmartDigraph &)
            NodeIt(const SmartDigraph &, const Node&)
            NodeIt& operator++()

        cppclass Arc:
            Arc()
            Arc(const Arc &)
            Arc(Invalid)
            bool operator==(Arc) const
            bool operator<=(Arc) const
            bool operator<(Arc) const

        cppclass ArcIt(Arc):
            ArcIt()
            ArcIt(Invalid)
            ArcIt(const SmartDigraph &)
            Arcit(const SmartDigraph &, Node &)
            ArcIt& operator++()

        cppclass OutArcIt(Arc):
            OutArcIt()
            OutArcIt(Invalid)
            OutArcIt(const SmartDigraph &)
            OutArcit(const SmartDigraph &, Node &)
            OutArcIt& operator++()

        cppclass InArcIt(Arc):
            InArcIt()
            InArcIt(Invalid)
            InArcIt(const SmartDigraph &)
            InArcit(const SmartDigraph &, Node &)
            InArcIt& operator++()

        Node baseNode(const OutArcIt &) const
        Node runningNode(const OutArcIt &) const
        Node baseNode(const InArcIt &) const
        Node runningNode(const InArcIt &) const
        
        Node addNode()
        Arc addArc(Node, Node)
        bool valid(Node) const
        bool valid(Arc) const
        Node split(Node)
        Node split(Node, bool)
        void clear()
        void reserveNode(int)
        void reserveArc(int)
        Node source(Arc) const
        Node target(Arc) const
        Node oppositeNode(const Node &, const Arc &) const
        int get_id "id"(Node) const
        int get_id "id"(Arc) const
        Node nodeFromId(int) const
        Arc arcFromId(int) const



cdef extern from "lemon/maps.h" namespace "lemon":
    cdef cppclass CrossRefMap[G, K, V]:
        CrossRefMap(const G&)
        void set(const K&, const V&)
        const V operator[](const K&) const
        K get "operator()"(const V&) const

