# boink/parsing.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

from cython.operator cimport dereference as deref
from libcpp.memory cimport make_unique

from boink.utils cimport _bstring

from khmer._oxli.parsing cimport Sequence

cdef class SplitPairedReader:

    def __init__(self, str left_filename, str right_filename,
                       int min_length=0, bool force_name_match=False):
        self._this = make_unique[_SplitPairedReader](_bstring(left_filename),
                                                     _bstring(right_filename),
                                                     min_length,
                                                     force_name_match)

    def __iter__(self):
        cdef _SequenceBundle bundle
        cdef object read_num = 0

        while not deref(self._this).is_complete():
            bundle = deref(self._this).next()
            yield (read_num, True,
                   Sequence._wrap(bundle.left),
                   Sequence._wrap(bundle.right))

            read_num += (<int>bundle.has_left + <int>bundle.has_right)
