from boink.pythonizors.utils import is_template_inst


def pythonize_boink(klass, name):
    is_dbg, _ = is_template_inst(name, 'dBG')
    if is_dbg:

        def add(self, item):
            if not isinstance(item, int) and len(item) < self.K:
                raise ValueError()
            elif isinstance(item, int) or len(item) == self.K:
                return self.insert(item)
            else:
                return self.insert_sequence(item)

        def wrap_query(func):
            def _get(self, item):
                return func(self, item)
            return _get

        def shallow_clone(self):
            return self.clone()

        def hashes(self, sequence):
            it = self.get_hash_iter(sequence)
            while not it.done():
                h = it.next()
                yield h

        def left_degree(self, kmer):
            self.set_cursor(kmer)
            return self.in_degree()

        def right_degree(self, kmer):
            self.set_cursor(kmer)
            return self.out_degree()

        def wrap_get(func):
            def _get(self):
                return func(self)
            return _get

        klass.add = add
        klass.get_hash = wrap_get(klass.get)
        klass.get = wrap_query(klass.query)
        klass.left_degree = left_degree
        klass.right_degree = right_degree
        klass.shallow_clone = shallow_clone
        klass.hashes = hashes 
        klass.K = property(klass.K)


