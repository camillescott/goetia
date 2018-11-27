# prometheus.pxd
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from libcpp.memory cimport weak_ptr, shared_ptr
from libcpp.string cimport string

cdef extern from "prometheus/collectable.h" namespace "prometheus" nogil:
    cdef cppclass _Collectable "prometheus::Collectable":
        pass

cdef extern from "prometheus/registry.h" namespace "prometheus" nogil:
    cdef cppclass _Registry "prometheus::Registry" (_Collectable):
        pass

cdef extern from "prometheus/exposer.h" namespace "prometheus" nogil:
    cdef cppclass _Exposer "prometheus::Exposer":
        _Exposer(const string&,
                 const string&)
        _Exposer(const string&)

        void RegisterCollectable(const weak_ptr[_Collectable]&)


cdef class Instrumentation:
    cdef shared_ptr[_Exposer] exposer;
    cdef shared_ptr[_Registry] registry;
