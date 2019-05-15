from boink import libboink
from boink.hashing import types as hasher_types
from boink.storage import get_storage_args, process_storage_args


def get_graph_args(parser):
    group = parser.add_argument_group('dBG')
    group.add_argument('-K', '--ksize', type=int, default=31)
    group.add_argument('--hasher', choices=list(hasher_types.keys()),
                       default='RollingHashShifter')

    get_storage_args(parser)

    return group


def process_graph_args(args):
    process_storage_args(args)

    args.graph_t = libboink.dBG[args.storage, getattr(libboink.hashing, args.hasher)]


def print_dBG_args(args):
    print('* dBG will be order', args.ksize, file=sys.stderr)
    print('* dBG will have underlying storage of', args.storage, file=sys.stderr)
    print('*','*' * 10, '*', sep='\n', file=sys.stderr)
