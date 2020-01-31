from boink.pythonizors.utils import is_template_inst

from cppyy import gbl

def pythonize_boink(klass, name):

    is_inst, _ = is_template_inst(name, 'FileProcessor')
    if is_inst:
        klass.advance.__release_gil__ = True
        klass.process.__release_gil__ = True

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

        def wrap_build(build_func):
            def wrapped(*args, fine_interval=10000,
                               medium_interval=100000,
                               coarse_interval=1000000,
                               verbose=False):
                return build_func(*args, fine_interval,
                                         medium_interval,
                                         coarse_interval,
                                         verbose)
            return wrapped

        
        klass.build = wrap_build(klass.build)
        klass.chunked_process = chunked_process

    if 'interval_state' in name:
        klass.__str__ = lambda self: '<interval_state fine={0} medium={1} coarse={2} end={3}'.format(self.fine, self.medium, self.coarse, self.end)
        def get(self):
            if self.fine:
                return 'fine'
            elif self.medium:
                return 'medium'
            elif self.coarse:
                return 'coarse'
            else:
                return 'end'
        klass.get = get

