from boink import boink as libboink

import argparse
import sys

types = [(libboink.storage.SparseppSetStorage, ()),
         (libboink.storage.BitStorage, (100000, 4)),
         (libboink.storage.ByteStorage, (100000, 4)),
         (libboink.storage.NibbleStorage, (100000, 4))]


def add_storage_args(parser):

    def storage_type(type_name):
        return getattr(libboink.storage, type_name)

    type_dict = {storage_t.__name__ : defaults for storage_t, defaults in types.storage_types}

    group = parser.add_argument_group('storage')

    group.add_argument('--storage', choices=list(type_dict.keys()), 
                        type=storage_type, default='SparseppSetParser')
    group.add_argument('-N', '--n_tables', default=4, type=int)
    group.add_argument('-x', '--max-tablesize', default=1e6, type=int)

    return group
