from cython.operator cimport dereference as deref

from boink.utils cimport _bstring, _ustring

cdef void _test():
    cdef DefaultShifter * shifter = new DefaultShifter(4)
    print(deref(shifter).set_cursor(_bstring('AAAA')))
    print(deref(shifter).get_cursor())
    print(deref(shifter).shift_left(b'C'))
    print(deref(shifter).get_cursor())

    print(deref(shifter).hash(_bstring('AAAA')))
    print(deref(shifter).hash(_bstring('CAAA')))

def test():
    _test()
