from boink.pythonizors.utils import *


def pythonize_boink_cdbg(klass, name):

    if name == 'compact_segment':

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

        klass.sequence = sequence
        klass.NULL = property(NULL)
        klass.start = property(start)
        klass.__repr__ = __repr__
        klass.__str__ = __str__

    cDBG_inst, template = is_template_inst(name, 'cDBG')
    if cDBG_inst:
        print()
        print('pythonizor:', name, klass)
        print('sub:', klass.Graph)
        klass.Graph.n_updates = property(klass.Graph.n_updates)
        klass.Graph.n_unodes = property(klass.Graph.n_unitig_nodes)
        klass.Graph.n_tags = property(klass.Graph.n_tags)
        klass.Graph.n_dnodes = property(klass.Graph.n_decision_nodes)
        klass.Graph.n_unitig_ends = property(klass.Graph.n_unitig_ends)
        print(dir(klass.Graph))

    compactor_inst, template = is_template_inst(name, 'StreamingCompactor')
    if compactor_inst:
        print()
        print('pythonizor: ', name, klass)
        klass.Compactor.Graph = property(lambda self: self.cdbg)


