import cppyy
from .initializor import initialize

__version__ = '0.15.4.1'

initialize('goetia', 'libgoetiaCppyy.so', 'goetia.map')
del initialize

from goetia import goetia as libgoetia
from cppyy import nullptr
