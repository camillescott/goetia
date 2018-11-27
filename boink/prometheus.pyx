# prometheus.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref
from libcpp.memory cimport make_shared

from boink.utils cimport _bstring

def is_port_in_use(port):
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

cdef class Instrumentation:

    def __cinit__(self, str port, str uri='/metrics'):
        _port = int(port)
        try:
            if is_port_in_use(_port):
                raise ValueError('Port {0} is in use.'.format(_port))
        except Exception as e:
            raise e

        _address = _bstring('127.0.0.1:' + port)
        _uri     = _bstring(uri)

        self.exposer  = make_shared[_Exposer](_address, _uri);
        self.registry = make_shared[_Registry]()
        cdef weak_ptr[_Registry] ref = weak_ptr[_Registry](self.registry)
        deref(self.exposer).RegisterCollectable(<weak_ptr[_Collectable]>ref)
