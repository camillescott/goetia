#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : traversal.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.10.2019

from boink import libboink
from boink.utils import copy_attrs
from cppyy.gbl import std
from cppyy import gbl

Traverse = libboink.Traverse
STATES = libboink.TraversalState


class Assembler:
    
    def __init__(self, graph):
        self.graph = graph
        self.graphptr = self.graph.__smartptr__().get()
        self.traverser = libboink.Traverse[type(graph)].dBG(graph.get_hasher())
        copy_attrs(libboink.Traverse[type(graph)],
                   self,
                   ['hash_type', 'shift_type', 'kmer_type'])
    
    def __getattr__(self, arg):
        attr = getattr(self.traverser, arg)
        if callable(attr) and not attr.__name__.startswith('__'):
            def wrapper(*args, **kwargs):
                return attr(self.graphptr, *args, **kwargs)
            return wrapper
    
    def set_cursor(self, kmer):
        self.traverser.set_cursor(kmer)

    @property
    def cursor(self):
        return self.traverser.get_cursor()

    @cursor.setter
    def cursor(self, seed):
        self.traverser.set_cursor(seed)
            
    def clear_seen(self):
        self.traverser.clear_seen()
    
    def assemble_right(self, seed):
        path = libboink.Path()
        mask = std.set[self.hash_type]()
        end_state = self.traverser.traverse_right(self.graphptr, seed, path, mask)
        return self.traverser.to_string(path), end_state

    def assemble_left(self, seed):
        path = libboink.Path()
        mask = std.set[self.hash_type]()
        end_state = self.traverser.traverse_left(self.graphptr, seed, path, mask)
        return self.traverser.to_string(path), end_state
    
    def assemble(self, seed):
        path = libboink.Path()
        mask = std.set[self.hash_type]()
        left_end_state, right_end_state = self.traverser.traverse(self.graphptr, seed, path, mask)
        return self.traverser.to_string(path), left_end_state, right_end_state
