import os

import cppyy
from .initializor import initialize

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'VERSION')) as fp:
    __version__ = fp.read().strip()

initialize('goetia', 'libgoetiaCppyy.so', 'goetia.map')
del initialize

from goetia import goetia as libgoetia
from cppyy import nullptr
