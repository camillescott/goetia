# cython: c_string_type=unicode, c_string_encoding=utf8, embedsignature=True

from cython cimport numeric

cdef class Weighted:

    def __cinit__(self, object item, float weight):
        self.item = item
        self.weight = weight

    def __richcmp__(Weighted x, Weighted y, int op):
        if op == 0:
            return x.weight < y.weight
        elif op == 1:
            return x.weight <= y.weight
        elif op == 2:
            return x.weight == y.weight
        elif op == 3:
            return x.weight != y.weight
        elif op == 4:
            return x.weight > y.weight
        else:
            return x.weight >= y.weight


cdef class MinHeap:

    def __cinit__(self, int capacity):
        self.data = [Weighted(None, 0) for _ in range(capacity)]
        self.n_items = 0
        self.capacity = capacity

    def values(self):
        for item in self.data:
            yield item.weight

    def keys(self):
        for item in self.data:
            yield item.item

    def iteritems(self):
        for item in self.data:
            yield item.item, item.weight

    cdef int left(self, int i):
        return 2 * i + 1

    cdef int right(self, int i):
        return 2 * i + 2

    cdef int parent(self, int i):
        return (i - 1) / 2 

    cdef Weighted root(self):
        return self.data[0]

    cpdef float min(self):
        if self.n_items < self.capacity:
            return min(self.data).weight
        return self.root().weight

    cdef void heapify(self, int i=-1):
        cdef int l, r
        cdef int m = i
        if i < 0:
            i = (self.capacity - 1) / 2
            while i >= 0:
                self.heapify(i)
                i = i - 1
        else:
            l = self.left(i)
            r = self.right(i)
            if l < self.capacity and (self.data[l] < self.data[i]):
                m = l
            if r < self.capacity and (self.data[r] < self.data[m]):
                m = r
            if m != i:
                self.data[i], self.data[m] = self.data[m], self.data[i]
                self.heapify(m)

    cpdef void insert(self, Weighted item):
        if self.n_items < self.capacity:
            self.data[self.n_items] = item
            self.n_items += 1
            if self.n_items == self.capacity:
                self.heapify()
        else:
            if item.weight > self.min():
                self.data[0] = item
                self.heapify(0)

        
