#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : pythonize_storage.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

from goetia.pythonizors.utils import make_tuple_from_type


def pythonize_goetia(klass, name):
    from goetia.utils import pretty_repr

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
                built = None
                if len(args) == 0:
                    built = build_func(klass.default_params)
                elif len(args) == 1 and isinstance(args[0], klass.Traits.params_type):
                    built = build_func(*args)
                else:
                    built = build_func(klass.make_params(*args))
                built.build_params = args
                return built
            return wrapped
        klass.build = wrap_build(klass.build)

        def describe(self):
            desc = f'{klass.__name__}\n'\
                   f'- Params: {self.build_params}\n'\
                   f'- Probabilistic: {bool(self.is_probabilistic)}\n'\
                   f'- Counting:      {bool(self.is_counting)}'
            return desc
        klass.describe = describe
