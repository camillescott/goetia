from boink.pythonizors.utils import is_template_inst
from boink.metadata import DATA_DIR


def pythonize_boink_hashing(klass, name):
    if name == 'UKHS':

        def get_kmers(W, K):
            import gzip

            valid_W = list(range(20, 210, 10))
            valid_K = list(range(7, 11))
            W = W - (W % 10)

            if not W in valid_W:
                raise ValueError('Invalid UKHS window size.')
            if not K in valid_K:
                raise ValueError('Invalid UKHS K.')

            filename = os.path.join(DATA_DIR,
                                    'res_{0}_{1}_4_0.txt.gz'.format(K, W))
            kmers = []
            with gzip.open(filename) as fp:
                for line in fp:
                    kmers.append(_bstring(line.strip()))

            return kmers 

        klass.get_kmers = staticmethod(get_kmers)
    
        
        def __uk_repr__(self):
            return '<Unikmer h={0} p={1}>'.format(self.hash, self.partition)

        def __uk_str__(self):
            return repr(self)

        klass.Unikmer.__repr__ = __uk_repr__
        klass.Unikmer.__str__  = __uk_str__


        def __binned_repr__(self):
            return '<BinnedKmer h={0} bin={1}>'.format(self.hash, self.unikmer)

        klass.BinnedKmer.__repr__ = __binned_repr__
        klass.BinnedKmer.__str__  = __binned_repr__
