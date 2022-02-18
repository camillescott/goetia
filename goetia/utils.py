#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : utils.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.10.2019

import collections
from dataclasses import is_dataclass
import itertools
import os
import re
import time
from typing import TypeVar, Type, Callable, List, Dict, Any


primitives = {
    'unsigned long': 'uint64_t',
    'long': 'uint64_t',
    'unsigned short': 'uint16_t'
}


def is_iterable(obj):
    return (
        isinstance(obj, collections.Iterable)
        and not isinstance(obj, str)
    )


def time_iterable(iterable):
    start = time.perf_counter()
    for result in iterable:
        elapsed = time.perf_counter() - start
        yield result, start, elapsed
        start = time.perf_counter()


class Counter:
    def __init__(self, val=0):
        self.count = val

    def __int__(self):
        return self.count

    def __add__(self, other):
        return Counter(self.count + other)

    def __iadd__(self, other):
        self.count += other
        return self

    def __eq__(self, other):
        return self.count == other

    def __ge__(self, other):
        return self.count >= other


def grouper(n, iterable):
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])


def find_common_basename(a, b):
    # it ain't elegant but it works
    a = os.path.basename(a)
    b = os.path.basename(b)
    for i in range(len(b), 0, -1):
        if a.find(b[:i]) != -1:
            return a[:i].strip('._-')
    return ''


def remove_fx_suffix(fn):
    return re.sub(r'(?:fq|FQ|fastq|FASTQ|fa|FA|fasta|FASTA)$', '', fn).strip('._-')


def check_trait(trait, type):
    return trait[type].value


def copy_attrs(src, dst, attrs):
    for attr in attrs:
        setattr(dst, attr, getattr(src, attr))


def set_typedef_attrs(dst_klass, attrs):
    for attr in attrs:
        setattr(dst_klass, attr, property(lambda self: getattr(type(self), attr)))


def pretty_repr(tname):
    namespaces = ['goetia', 'hashing', 'sequences', 'kmers', 'storage', 'cdbg', 'parsing', 'reporting', 'events', 'metrics']
    pretty = tname.__name__.replace(' >', '>').replace('  >', '>')
    for namespace in namespaces:
        pretty = pretty.replace(namespace + '::', '')
    for old, new in primitives.items():
        pretty = pretty.replace(old, new)
    pretty = '(' + pretty + ')'
    return pretty


class Namespace(): pass


#
# Dataclass utilities taken from: https://gist.github.com/mikeholler/4be180627d3f8fceb55704b729464adb
#

_T = TypeVar("_T")
_Self = TypeVar("_Self")
_VarArgs = List[Any]
_KWArgs = Dict[str, Any]

def _kwarg_only_init_wrapper(
        self: _Self,
        init: Callable[..., None],
        *args: _VarArgs,
        **kwargs: _KWArgs
) -> None:
    if len(args) > 0:
        raise TypeError(
            f"{type(self).__name__}.__init__(self, ...) only allows keyword arguments. Found the "
            f"following positional arguments: {args}"
        )
    init(self, **kwargs)


def _positional_arg_only_init_wrapper(
        self: _Self,
        init: Callable[..., None],
        *args: _VarArgs,
        **kwargs: _KWArgs
) -> None:
    if len(kwargs) > 0:
        raise TypeError(
            f"{type(self).__name__}.__init__(self, ...) only allows positional arguments. Found "
            f"the following keyword arguments: {kwargs}"
        )
    init(self, *args)


def require_kwargs(cls: Type[_T]) -> Type[_T]:
    """
    Force a dataclass's init function to only work if called with keyword arguments.
    If parameters are not positional-only, a TypeError is thrown with a helpful message.
    This function may only be used on dataclasses.
    This works by wrapping the __init__ function and dynamically replacing it. Therefore,
    stacktraces for calls to the new __init__ might look a bit strange. Fear not though,
    all is well.
    Note: although this may be used as a decorator, this is not advised as IDEs will no longer
    suggest parameters in the constructor. Instead, this is the recommended usage::
        from dataclasses import dataclass
        @dataclass
        class Foo:
            bar: str
        require_kwargs_on_init(Foo)
    """

    if cls is None:
        raise TypeError("Cannot call with cls=None")
    if not is_dataclass(cls):
        raise TypeError(
            f"This decorator only works on dataclasses. {cls.__name__} is not a dataclass."
        )

    original_init = cls.__init__

    def new_init(self: _Self, *args: _VarArgs, **kwargs: _KWArgs) -> None:
        _kwarg_only_init_wrapper(self, original_init, *args, **kwargs)

    # noinspection PyTypeHints
    cls.__init__ = new_init  # type: ignore

    return cls


def require_positional_args(cls: Type[_T]) -> Type[_T]:
    """
    Force a dataclass's init function to only work if called with positional arguments.
    If parameters are not positional-only, a TypeError is thrown with a helpful message.
    This function may only be used on dataclasses.
    This works by wrapping the __init__ function and dynamically replacing it. Therefore,
    stacktraces for calls to the new __init__ might look a bit strange. Fear not though,
    all is well.
    Note: although this may be used as a decorator, this is not advised as IDEs will no longer
    suggest parameters in the constructor. Instead, this is the recommended usage::
        from dataclasses import dataclass
        @dataclass
        class Foo:
            bar: str
        require_positional_args_on_init(Foo)
    """

    if cls is None:
        raise TypeError("Cannot call with cls=None")
    if not is_dataclass(cls):
        raise TypeError(
            f"This decorator only works on dataclasses. {cls.__name__} is not a dataclass."
        )

    original_init = cls.__init__

    def new_init(self: _Self, *args: _VarArgs, **kwargs: _KWArgs) -> None:
        _positional_arg_only_init_wrapper(self, original_init, *args, **kwargs)

    # noinspection PyTypeHints
    cls.__init__ = new_init  # type: ignore

    return cls
