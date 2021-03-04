#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : pythonize_storage.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

from goetia.pythonizors.utils import make_tuple_from_type


def pythonize_goetia_storage(klass, name):

    if name.strip(' ').endswith('Storage'):

        klass.is_counting = bool(klass.Traits.is_counting)
        klass.is_probabilistic = bool(klass.Traits.is_probabilistic)
        klass.params_type = klass.Traits.params_type
        klass.default_params = klass.Traits.default_params

        def make_params(*args):
            return make_tuple_from_type(klass.params_type, *args)
        klass.make_params = staticmethod(make_params)

        def wrap_build(build_func):
            def wrapped(*args):
                if len(args) == 0:
                    return build_func(klass.default_params)
                if len(args) == 1 and isinstance(args[0], klass.Traits.params_type):
                    return build_func(*args)
                return build_func(klass.make_params(*args))
            return wrapped
        klass.build = wrap_build(klass.build)