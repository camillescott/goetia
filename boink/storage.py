from boink import libboink
from boink.utils import check_trait

import argparse
import sys

_types = [(libboink.storage.SparseppSetStorage, tuple()),
          (libboink.storage.BitStorage, (100000, 4)),
          (libboink.storage.ByteStorage, (100000, 4)),
          (libboink.storage.NibbleStorage, (100000, 4))]

types = {storage_t.__name__ : defaults for storage_t, defaults in _types}


SparseppSetStorage = libboink.storage.SparseppSetStorage
BitStorage         = libboink.storage.BitStorage
ByteStorage        = libboink.storage.ByteStorage
NibbleStorage      = libboink.storage.NibbleStorage

def is_counting(klass):
    return check_trait(libboink.storage.is_counting, klass)

def is_probabilistic(klass):
    return check_trait(libboink.storage.is_probabilistic, klass)


def get_storage_args(parser):
    if 'storage' in [g.title for g in parser._action_groups]:
        return None

    group = parser.add_argument_group('storage')

    group.add_argument('--storage', choices=list(types.keys()), 
                       default='SparseppSetStorage')
    group.add_argument('-N', '--n_tables', default=4, type=int)
    group.add_argument('-x', '--max-tablesize', default=int(1e5), type=int)

    return group


def process_storage_args(args):
    if hasattr(args, 'storage_args'):
        return

    args.storage = getattr(libboink.storage, args.storage)

    args.storage_args = ()
    if args.storage is not libboink.storage.SparseppSetStorage:
        args.storage_args = (args.max_tablesize, args.n_tables)

