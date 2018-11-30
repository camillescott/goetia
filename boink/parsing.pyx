# boink/parsing.pyx
# Copyright (C) 2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.

import itertools

cimport cython
from cython.operator cimport dereference as deref
from libcpp.memory cimport make_unique

from boink.utils cimport _bstring

PAIRING_MODES = ('split', 'interleaved', 'single')


def grouper(n, iterable):
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])


@cython.freelist(100)
cdef class Sequence:

    def __cinit__(self, name=None, sequence=None,
                        quality=None, description=None,
                        cleaned_seq=None):

        if name is not None and sequence is not None:
            self._obj.sequence = _bstring(sequence)
            self._obj.name = _bstring(name)
            if description is not None:
                self._obj.description = _bstring(description)
            if quality is not None:
                self._obj.quality = _bstring(quality)
            if cleaned_seq is not None:
                self._obj.cleaned_seq = _bstring(cleaned_seq)
            else:
                self._obj.cleaned_seq = self._obj.sequence

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'Sequence(name="{0}", sequence="{1}")'.format(self.name, self.sequence)

    def __len__(self):
        return self._obj.sequence.length()

    def __richcmp__(x, y, op):
        if op == 2:
            return x.name == y.name and x.sequence == y.sequence
        else:
            raise NotImplementedError('Operator not available')

    def kmers(self, int K):
        cdef int i = 0
        cdef unicode sequence = self.sequence
        for i in range(0, len(self)-K+1):
            yield sequence[i:i+K]

    def __getitem__(self, x):
        # Definitely optimize this.
        return self.sequence[x]

    @property
    def name(self):
        cdef unicode name = self._obj.name
        return self._obj.name if name else None

    @property
    def sequence(self):
        cdef unicode sequence = self._obj.sequence
        return self._obj.sequence if sequence else None

    @property
    def description(self):
        cdef unicode description = self._obj.description
        return description if description else None

    @property
    def quality(self):
        cdef unicode quality = self._obj.quality
        return quality if quality else None

    @property
    def cleaned_seq(self):
        cdef unicode cleaned_seq = self._obj.cleaned_seq
        return cleaned_seq if cleaned_seq else None

    @staticmethod
    def from_screed_record(record):
        cdef Sequence seq = Sequence(name=record.name,
                                     sequence=record.sequence)
        if hasattr(record, 'quality'):
            seq._obj.quality = _bstring(record.quality)

        for attr in ('annotations', 'description'):
            if hasattr(record, attr):
                seq._obj.description = _bstring(getattr(record, attr))

        return seq

    @staticmethod
    cdef Sequence _wrap(_Sequence cseq):
        cdef Sequence seq = Sequence()
        seq._obj = cseq
        return seq


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
