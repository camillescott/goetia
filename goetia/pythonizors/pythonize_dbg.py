from goetia.pythonizors.utils import is_template_inst


def pythonize_goetia(klass, name):
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

        def wrap_vector_ret(func):
            def wrapped(self, *args, **kwargs):
                return [item for item in func(self, *args, **kwargs)]
            return wrapped

        def wrap_walk(func):
            pass

        klass.add = add
        klass.get_hash = wrap_get(klass.get)
        klass.get = wrap_query(klass.query)
        klass.left_degree = left_degree
        klass.right_degree = right_degree
        klass.shallow_clone = shallow_clone
        klass.hashes = hashes

        klass.left_extensions = wrap_vector_ret(klass.left_extensions)
        klass.right_extensions = wrap_vector_ret(klass.right_extensions)
        #klass.filter_nodes = wrap_vector_ret(klass.filter_nodes)
        #klass.in_neighbors = wrap_vector_ret(klass.in_neighbors)
        #klass.out_neighbors = wrap_vector_ret(klass.out_neighbors)

    is_pdbg, _ = is_template_inst(name, 'PdBG')
    if is_pdbg:
        def wrap_build(build_func):
            def wrapped(W, K, *storage_args, ukhs=None):
                if ukhs is None:
                    ukhs = klass.ukhs_type.load(W, K)
                if len(storage_args) == 0:
                    params = klass.base_storage_type.default_params
                elif len(storage_args) == 1 and \
                     isinstance(storage_args[0], klass.base_storage_type.params_type):
                    params = storage_args[0]
                else:
                    params = klass.base_storage_type.make_params(*storage_args)
                return build_func(W, K, ukhs, params)
            return wrapped

        klass.build = wrap_build(klass.build)

