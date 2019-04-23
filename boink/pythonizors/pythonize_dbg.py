from boink.pythonizors.utils import is_template_inst


def pythonize_boink(klass, name):
    if is_template_inst(name, 'dBG'):

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
            lneighbors = self.left_neighbors(kmer)
            return len(lneighbors)

        def right_degree(self, kmer):
            rneighbors = self.right_neighbors(kmer)
            return len(rneighbors)

        klass.add = add
        klass.get = get
        klass.left_degree = left_degree
        klass.right_degree = right_degree
        klass.shallow_clone = shallow_clone
        klass.hashes = hashes
