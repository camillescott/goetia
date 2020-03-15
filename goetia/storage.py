from goetia import libgoetia
from goetia.utils import check_trait

import argparse
import sys

_types = [(libgoetia.storage.SparseppSetStorage, tuple()),
          (libgoetia.storage.BitStorage, (100000, 4)),
          (libgoetia.storage.ByteStorage, (100000, 4)),
          (libgoetia.storage.NibbleStorage, (100000, 4))]

types = {storage_t.__name__.replace(' ', '') : defaults for storage_t, defaults in _types}


SparseppSetStorage = libgoetia.storage.SparseppSetStorage
BitStorage         = libgoetia.storage.BitStorage
ByteStorage        = libgoetia.storage.ByteStorage
NibbleStorage      = libgoetia.storage.NibbleStorage

count_t            = libgoetia.storage.count_t


def is_counting(klass):
    return check_trait(libgoetia.storage.is_counting, klass)


def is_probabilistic(klass):
    return check_trait(libgoetia.storage.is_probabilistic, klass)


def get_storage_args(parser, default='SparseppSetStorage'):
    if 'storage' in [g.title for g in parser._action_groups]:
        return None

    group = parser.add_argument_group('storage')

    group.add_argument('--storage', choices=list(types.keys()), 
                       default=default)
    group.add_argument('-N', '--n_tables', default=4, type=int)
    group.add_argument('-x', '--max-tablesize', default=1e5, type=float)

    return group


def process_storage_args(args):
    if hasattr(args, 'storage_args'):
        return

    args.storage = getattr(libgoetia.storage, args.storage)

    args.storage_args = ()
    if args.storage is not libgoetia.storage.SparseppSetStorage:
        args.max_tablesize = int(args.max_tablesize)
        args.storage_args = (args.max_tablesize, args.n_tables)

