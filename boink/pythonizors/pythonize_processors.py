from boink.pythonizors.utils import is_template_inst

from cppyy import gbl

def pythonize_boink(klass, name):

    is_inst, _ = is_template_inst(name, 'FileProcessor')
    if is_inst:
        def chunked_process(self, file, right_file=None, alphabet=None):
            if type(file) in (str, bytes):
                from boink.parsing import FastxParser, SplitPairedReader
                from boink.alphabets import DNA_SIMPLE

                if alphabet is None:
                    alphabet = DNA_SIMPLE
                parser_type = FastxParser[alphabet]

                if right_file is None:
                    parser = parser_type.build(file)
                else:
                    parser = SplitPairedReader[parser_type].build(file,
                                                                  right_file)
            else:
                parser = file
            
            while True:
                state = self.advance(parser)
                yield self.n_reads(), parser.n_skipped(), state
                if state.end:
                    break

        klass.chunked_process = chunked_process

    if 'interval_state' in name:
        klass.__str__ = lambda self: '<interval_state fine={0} medium={1} coarse={2} end={3}'.format(self.fine, self.medium, self.coarse, self.end)

