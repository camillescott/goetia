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

        def get(self, item):
            return self.query(item)

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

        klass.add = add
        klass.get = get
        klass.left_degree = left_degree
        klass.right_degree = right_degree
        klass.shallow_clone = shallow_clone
        klass.hashes = hashes 
        klass.K = property(klass.K)


