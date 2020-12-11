from goetia import libgoetia
from goetia.utils import check_trait

typenames = [(t, t.__name__.replace(' ', '')) for t in [libgoetia.storage.SparseppSetStorage,
                                                        libgoetia.storage.BitStorage,
                                                        libgoetia.storage.ByteStorage,
                                                        libgoetia.storage.NibbleStorage]]

types = [_type for _type, _name in typenames]

for hasher_t, name in typenames:
    globals()[name] = hasher_t


count_t = libgoetia.storage.count_t
StorageTraits = libgoetia.storage.StorageTraits



def get_storage_args(parser, default='SparseppSetStorage'):
    if 'storage' in [g.title for g in parser._action_groups]:
        return None

    group = parser.add_argument_group('storage')

    group.add_argument('--storage', choices=[name for _, name in typenames], 
                       default=default)
    group.add_argument('-N', '--n_tables', default=4, type=int)
    group.add_argument('-x', '--max-tablesize', default=1e8, type=float)

    return group


def process_storage_args(args):
    if hasattr(args, 'storage_args'):
        return

    args.storage = getattr(libgoetia.storage, args.storage)

    args.storage_args = ()
    if args.storage is not libgoetia.storage.SparseppSetStorage:
        args.max_tablesize = int(args.max_tablesize)
        args.storage_args = (args.max_tablesize, args.n_tables)

