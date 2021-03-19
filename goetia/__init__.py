import os

import cppyy
from .initializor import initialize

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'VERSION')) as fp:
    __version__ = fp.read().strip()

initialize('goetia', 'libgoetiaCppyy.so', 'goetia.map')
del initialize

from goetia import goetia as libgoetia
from cppyy import nullptr


splash = f'''
                                   d8,          
                            d8P   `8P           
                         d888888P               
 d888b8b   d8888b  d8888b  ?88'    88b d888b8b  
d8P' ?88  d8P' ?88d8b_,dP  88P     88Pd8P' ?88  
88b  ,88b 88b  d8888b      88b    d88 88b  ,88b 
`?88P'`88b`?8888P'`?888P'  `?8b  d88' `?88P'`88b
       )88                                      
      ,88P                        {__version__}
  `?8888P                                       

'''
