#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : storage.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 12.03.2021

import math

from goetia import libgoetia
from goetia.utils import check_trait

typenames = [(t, t.__name__.replace(' ', '')) for t in [libgoetia.SparseppSetStorage,
                                                        libgoetia.PHMapStorage,
                                                        libgoetia.BitStorage,
                                                        libgoetia.ByteStorage,
                                                        libgoetia.NibbleStorage,
                                                        libgoetia.QFStorage,
                                                        libgoetia.BTreeStorage,
                                                        libgoetia.HLLStorage]]

types = [_type for _type, _name in typenames]

for hasher_t, name in typenames:
    globals()[name] = hasher_t


count_t = libgoetia.count_t
StorageTraits = libgoetia.StorageTraits


def get_storage_args(parser, default='SparseppSetStorage',
                     group_name='storage'):
    if 'storage' in [g.title for g in parser._action_groups]:
        return None

    group = parser.add_argument_group(group_name)

    group.add_argument('-S', '--storage',
                       choices=[name for _, name in typenames], 
                       default=default)
    group.add_argument('-N', '--n_tables',
                       default=4, type=int)
    group.add_argument('-x', '--max-tablesize',
                       default=1e8, type=float)
    group.add_argument('-E', '--error-rate',
                       default=0.02, type=float)

    return group


def process_storage_args(args):
    if hasattr(args, 'storage_args'):
        return

    args.storage = getattr(libgoetia, args.storage)

    args.storage_args = ()
    if args.storage in (libgoetia.BitStorage,
                        libgoetia.ByteStorage,
                        libgoetia.NibbleStorage):

        args.max_tablesize = int(args.max_tablesize)
        args.storage_args = (args.max_tablesize, args.n_tables)

    elif args.storage is libgoetia.HLLStorage:
        args.storage_args = (float(args.error_rate), )

    elif args.storage is libgoetia.QFStorage:
        args.storage_args = (int(math.ceil(math.log2(args.max_tablesize))), )

    else:
        args.storage_args = tuple()


def get_partitioned_storage_args(parser):
    group = get_storage_args(parser, 'partitioned storage')

