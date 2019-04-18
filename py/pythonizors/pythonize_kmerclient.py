from cppyy import gbl


def pythonize_boink_kmers(klass, name):

    if name == 'KmerClient':

        klass.K = property(klass.K)
