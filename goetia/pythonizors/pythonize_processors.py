from goetia.pythonizors.utils import is_template_inst

from cppyy import gbl

def pythonize_goetia(klass, name):

    is_inst, _ = is_template_inst(name, 'FileProcessor')
    if is_inst:
        klass.advance.__release_gil__ = True
        klass.process.__release_gil__ = True

        def chunked_process(self, file, right_file=None):
            if type(file) in (str, bytes):
                from goetia.parsing import FastxParser, SplitPairedReader

                parser_type = FastxParser[type(self).alphabet]

                if right_file is None:
                    parser = parser_type.build(file)
                else:
                    parser = SplitPairedReader[parser_type].build(file,
                                                                  right_file)
            else:
                parser = file
            
            advancing = True
            prev_n_sequences = 0
            while advancing:
                n_sequences, time_elapsed, advancing = self.advance(parser)
                if n_sequences > prev_n_sequences:
                    yield n_sequences, time_elapsed, parser.n_skipped()
                prev_n_sequences = n_sequences


        def wrap_build(build_func):
            def wrapped(*args, interval=10000,
                               verbose=False):
                return build_func(*args, interval,
                                         verbose)
            return wrapped

        
        klass.build = wrap_build(klass.build)
        klass.chunked_process = chunked_process
