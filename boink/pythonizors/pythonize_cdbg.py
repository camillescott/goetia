from boink.pythonizors.utils import *


def pythonize_boink_cdbg(klass, name):

    cDBG_inst, template = is_template_inst(name, 'cDBG')

    if cDBG_inst:

        #
        # CompactNode pythonizations
        #
    
        klass.CompactNode.__len__ = lambda self: self.length()
        klass.CompactNode.meta = property(klass.CompactNode.meta)

        #
        # UnitigNode pythonizations
        #
    
        klass.UnitigNode.left_end = property(klass.UnitigNode.left_end)
        klass.UnitigNode.right_end = property(klass.UnitigNode.right_end)

        #
        # Graph pythonizations
        #

        klass.Graph.n_updates = property(klass.Graph.n_updates)
        klass.Graph.n_unodes = property(klass.Graph.n_unitig_nodes)
        klass.Graph.n_tags = property(klass.Graph.n_tags)
        klass.Graph.n_dnodes = property(klass.Graph.n_decision_nodes)
        klass.Graph.n_unitig_ends = property(klass.Graph.n_unitig_ends)

        klass.Graph.query_cnode = convert_nullptr(klass.Graph.query_cnode)
        klass.Graph.query_dnode = convert_nullptr(klass.Graph.query_dnode)
        klass.Graph.query_unode_end = convert_nullptr(klass.Graph.query_unode_end)
        klass.Graph.query_unode_tag = convert_nullptr(klass.Graph.query_unode_tag)
        klass.Graph.query_unode_id = convert_nullptr(klass.Graph.query_unode_id)


    compactor_inst, template = is_template_inst(name, 'StreamingCompactor')
    if compactor_inst:
        klass.Compactor.Graph = property(lambda self: self.cdbg)

        #
        # compact_segment pythonizations
        #

        def sequence(self, sequence: str) -> str:
            """sequence

            Args:
                Reference sequence the segment is from.

            Returns:
                str: Substring of the segment.
            """

            if self.NULL:
                return ''
            return sequence[self.start_pos : self.start_pos + self.length]

        def NULL(self) -> bool:
            """NULL

            Returns:
                bool: Whether the segment is a null segment.
            """
            return self.is_null()

        def start(self) -> int:
            return self.start_pos

        def __repr__(self):
            return '<Segment start={} length={} left_anchor={} '\
                     'right_anchor={} is_decision_kmer={} is_null={}>'.format(self.start_pos,
                                                                              self.length,
                                                                              self.left_anchor,
                                                                              self.right_anchor,
                                                                              self.is_decision_kmer,
                                                                              self.is_null())

        def __str__(self):
            return repr(self)

        klass.compact_segment.sequence = sequence
        klass.compact_segment.NULL = property(NULL)
        klass.compact_segment.start = property(start)
        klass.compact_segment.__repr__ = __repr__
        klass.compact_segment.__str__ = __str__
