#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : utils.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.10.2019

import itertools
import os
import re

primitives = {
    'unsigned long': 'uint64_t',
    'long': 'uint64_t',
    'unsigned short': 'uint64_t'
}


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