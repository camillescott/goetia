# boink/utils.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cpython.version cimport PY_MAJOR_VERSION
from math import sqrt


cdef bytes _bstring(s):
    if not isinstance(s, (basestring, bytes)):
        raise TypeError("Requires a string-like sequence")

    if isinstance(s, unicode):
        s = s.encode('utf-8')
    return s


cdef unicode _ustring(s):
    if type(s) is unicode:
        # fast path for most common case(s)
        return <unicode>s
    elif PY_MAJOR_VERSION < 3 and isinstance(s, bytes):
        # only accept byte strings in Python 2.x, not in Py3
        return (<bytes>s).decode('UTF-8')
    elif isinstance(s, unicode):
        # an evil cast to <unicode> might work here in some(!) cases,
        # depending on what the further processing does.  to be safe,
        # we can always create a copy instead
        return unicode(s)
    else:
        raise TypeError(...)


cpdef bool is_str(object s):
    return isinstance(s, (basestring, bytes))


cpdef bool is_num(object n):
    return isinstance(n, (int, long))


cpdef bool is_prime(uint64_t n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    cdef uint64_t i
    for i in range(3, sqrt(n), 2):
        if n % i == 0:
            return False

    return True;


cdef vector[uint64_t] get_n_primes_near_x(uint32_t n, uint64_t x):
    cdef vector[uint64_t] primes
    if x == 1:
        primes.push_back(1);
        return primes

    cdef uint64_t i = x - 1
    if i % 2 == 0:
        i -= 1

    while primes.size() != n:
        if is_prime(i):
            primes.push_back(i)

        if i == 1:
            break

        i -= 2

    return primes

