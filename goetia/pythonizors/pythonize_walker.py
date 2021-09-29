#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : pythonize_walker.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 22.01.2020


from cppyy import gbl
from cppyy.gbl import std

from goetia.pythonizors.utils import is_template_inst


def pythonize_goetia(klass, name):
    is_walker, _ = is_template_inst(name, 'UnitigWalker')
    if is_walker:

        # def wrap_walk(wrapped):
        #     def _walk(self, seed=None, mask=None):
        #         if mask is None:
        #             mask = std.set[type(self).hash_type]()
        #         if seed is None:
        #             return wrapped(self, mask)
        #         else:
        #             return wrapped(self, seed, mask)
        #     return _walk

        # klass.walk_left = wrap_walk(klass.walk_left)
        # klass.walk_right = wrap_walk(klass.walk_right)
        # klass.walk = wrap_walk(klass.walk)

        klass.cursor = property(klass.get_cursor)
        klass.cursor = klass.cursor.setter(lambda self, seed: self.set_cursor(seed))

        klass.Walk[gbl.goetia.hashing.DIR_LEFT].__str__ = klass.Walk[gbl.goetia.hashing.DIR_LEFT].to_string
        klass.Walk[gbl.goetia.hashing.DIR_RIGHT].__str__ = klass.Walk[gbl.goetia.hashing.DIR_RIGHT].to_string
