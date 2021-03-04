from goetia.pythonizors.utils import is_template_inst
from goetia.utils import set_typedef_attrs
from goetia.metadata import DATA_DIR
from cppyy.gbl import std

UKHS_CACHE = dict()

def pythonize_goetia_hashing(klass, name):

    ukhs_inst, template = is_template_inst(name, 'UKHS')
    if ukhs_inst:
        
        def parse_unikmers(W, K):
            import gzip
            import os

            valid_W = list(range(20, 210, 10))
            valid_K = list(range(7, 11))
            W = W - (W % 10)

            if not W in valid_W:
                raise ValueError('Invalid UKHS window size.')
            if not K in valid_K:
                raise ValueError('Invalid UKHS K.')

            filename = os.path.join(DATA_DIR,
                                    'res_{0}_{1}_4_0.txt.gz'.format(K, W))
            unikmers = std.vector[std.string]()
            with gzip.open(filename, 'rt') as fp:
                for line in fp:
                    unikmers.push_back(line.strip())

            return unikmers

        klass.parse_unikmers = staticmethod(parse_unikmers)

        def load(W, K):
            key = (W, K, klass.__name__)
            if key in UKHS_CACHE:
                return UKHS_CACHE[key]
            else:
                unikmers = parse_unikmers(W, K)
                ukhs = klass.build(W, K, unikmers)
                UKHS_CACHE[key] = ukhs
                return ukhs

        klass.load = staticmethod(load)

    shifter_inst, template = is_template_inst(name, 'HashShifter')
    if shifter_inst:
        def __getattr__(self, arg):
            attr = getattr(type(self), arg)
            if not attr.__name__.startswith('__'):
                return attr

        klass.__getattr__ = __getattr__

        if 'Unikmer' in name:
            def build(W, K):
                ukhs_type = klass.ukhs_type
                ukhs = ukhs_type.load(W, K)
                shifter = klass(W, K, ukhs.__smartptr__())
                return shifter

            klass.build = staticmethod(build)


        #set_typedef_attrs(klass, ['alphabet', 'hash_type', 'value_type', 'kmer_type'])

    for check_name in ['Hash', 'Canonical', 'Wmer',
                       'Kmer', 'Shift', 'Partitioned']:

        is_inst, _ = is_template_inst(name, check_name)
        if is_inst:

            klass.value = property(klass.value)
            klass.__lt__ = lambda self, other: self.value < other.value
            klass.__le__ = lambda self, other: self.value <= other.value
            klass.__gt__ = lambda self, other: self.value > other.value
            klass.__ge__ = lambda self, other: self.value >= other.value
            klass.__ne__ = lambda self, other: self.value != other.value
            klass.__repr__ = klass.__str__
            klass.__hash__ = lambda self: hash(self.value)
