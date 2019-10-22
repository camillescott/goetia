from boink.pythonizors.utils import is_template_inst

from cppyy import gbl

def pythonize_boink(klass, name):

    if is_template_inst('FileProcessor', name):
        def chunked_process(self, filename, right_filename=None):
            from boink.parsing import FastxReader, SplitPairedReader
            if right_filename is None:
                parser = FastxReader.build(filename)
            else:
                parser = SplitPairedReader(filename,
                                           right_filename)
            while True:
                state = self.advance(parser)
                yield self.n_reads(), state
                if state.end:
                    break

        klass.chunked_process = chunked_process

    if 'interval_state' in name:
        klass.__str__ = lambda self: '<interval_state fine={0} medium={1} coarse={2} end={3}'.format(self.fine, self.medium, self.coarse, self.end)

