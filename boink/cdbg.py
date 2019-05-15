from boink import libboink
from boink.cli import CommandRunner, get_output_interval_args
from boink.dbg import get_graph_args, process_graph_args
from boink.parsing import get_pairing_args, iter_fastx_inputs
from boink.prometheus import Instrumentation, get_prometheus_args
from boink.serialization import cDBGSerialization
from boink.metadata import CUR_TIME

import os


def get_cdbg_args(parser):
    default_prefix = 'boink.build-cdbg.' + CUR_TIME
    parser.default_prefix = default_prefix
    group = parser.add_argument_group('cDBG')

    group.add_argument('--results-dir',
                        default=default_prefix)

    group.add_argument('--normalize',
                        type=int,
                        nargs='?',
                        const=10)

    group.add_argument('--save-cdbg',
                        metavar='PREFIX.<format>',
                        nargs='?',
                        const='boink.cdbg.graph')
    group.add_argument('--save-cdbg-format',
                        nargs='+',
                        choices=cDBGSerialization.FORMATS,
                        default=['gfa1'])

    group.add_argument('--track-cdbg-stats',
                        metavar='FILE_NAME.csv',
                        nargs='?',
                        const='boink.cdbg.stats.csv')

    group.add_argument('--track-cdbg-history',
                        metavar='FILENAME.graphml',
                        nargs='?',
                        const='boink.cdbg.history.graphml')

    group.add_argument('--track-cdbg-components',
                        metavar='FILE_NAME.csv',
                        nargs='?',
                        const='boink.cdbg.components.csv')
    group.add_argument('--component-sample-size',
                        type=int,
                        default=10000)

    group.add_argument('--track-cdbg-unitig-bp',
                        metavar='FILENAME.csv',
                        nargs='?',
                        const='boink.cdbg.unitigs.bp.csv')

    group.add_argument('--unitig-bp-bins',
                        nargs='+',
                        type=int)

    group.add_argument('--validate',
                        metavar='FILENAME.csv',
                        nargs='?',
                        const='boink.cdbg.validation.csv')

    return group


def process_cdbg_args(args):
    def join(p):
        return p if p is None else os.path.join(args.results_dir, p)
    args.track_cdbg_stats =      join(args.track_cdbg_stats)
    args.track_cdbg_components = join(args.track_cdbg_components)
    args.track_cdbg_history =    join(args.track_cdbg_history)
    args.save_cdbg =             join(args.save_cdbg)
    args.track_cdbg_unitig_bp =  join(args.track_cdbg_unitig_bp)




def print_cdbg_args(args):
    print('* cDBG Params', file=sys.stderr)
    print('* Directory: ', args.results_dir, file=sys.stderr)
    if args.save_cdbg:
        print('* Saving cDBG every {0} sequences with file prefix {1}'.format(args.coarse_interval,
                                                                              args.save_cdbg),
              file=sys.stderr)
        print('* cDBG save formats: {0}'.format(', '.join(args.save_cdbg_format)))
    if args.track_cdbg_stats:
        print('* Tracking cDBG stats and reporting every {0} sequences'.format(args.fine_interval),
              file=sys.stderr)
        print('* Saving tracking information to', args.track_cdbg_stats, file=sys.stderr)
    if args.track_cdbg_history:
        print('* Tracking cDBG history and saving to', args.track_cdbg_history, file=sys.stderr)
    if args.validate:
        print('* cDBG will be validated on completion and results saved to', args.validate,
              file=sys.stderr)
    print('*', '*' * 10, '*', sep='\n', file=sys.stderr)


class cDBGRunner(CommandRunner):

    def __init__(self, parser):
        get_graph_args(parser)
        get_cdbg_args(parser)
        get_prometheus_args(parser)
        get_output_interval_args(parser)
        group = get_pairing_args(parser)
        group.add_argument('-o', dest='output_filename', default='/dev/stdout')
        group.add_argument('-i', dest='inputs', nargs='+', required=True)

        super().__init__(parser)

    def postprocess_args(self, args):
        process_graph_args(args)
        process_cdbg_args(args)

    def setup(self, args):
        os.makedirs(args.results_dir, exist_ok=True)

        self.instrumentation = Instrumentation(args.port,
                                               expose=(args.port != None))

        self.dbg_t       = args.graph_t
        self.dbg         = args.graph_t.build(args.ksize,
                                              *args.storage_args)

        self.cdbg_t      = libboink.cdbg.cDBG[type(self.dbg)]

        self.compactor_t = libboink.cdbg.StreamingCompactor[type(self.dbg)]
        self.compactor = self.compactor_t.Compactor.build(self.dbg,
                                                          self.instrumentation.Registry)

        self.processor = self.compactor_t.Processor[libboink.parsing.FastxReader].build(self.compactor,
                                                                                        args.fine_interval,
                                                                                        args.medium_interval,
                                                                                        args.coarse_interval)
        if args.track_cdbg_stats:
            self.cdbg_reporter = self.compactor_t.Reporter.build(args.track_cdbg_stats,
                                                                 self.compactor)
            self.processor.register_listener(self.cdbg_reporter)

        if args.track_cdbg_unitig_bp:
            if args.unitig_bp_bins is None:
                bins = [args.ksize, 100, 200, 500, 1000]
            else:
                bins = args.unitig_bp_bins
            self.unitig_reporter = self.cdbg_t.UnitigReporter.build(args.track_cdbg_unitig_bp,
                                                                    self.compactor.cdbg.__smartptr__(),
                                                                    bins)
            self.processor.register_listener(self.unitig_reporter)

        if args.track_cdbg_components:
            self.components = self.cdbg_t.ComponentReporter.build(args.track_cdbg_components,
                                                                  self.compactor.cdbg.__smartptr__(),
                                                                  args.component_sample_size,
                                                                  self.instrumention.Registry)
            self.processor.register_listener(self.components)

        self.writers = []
        if args.save_cdbg:
            for cdbg_format in args.save_cdbg_format:
                writer = self.cdbg_t.Writer.build(args.save_cdbg + '.' + cdbg_format,
                                                  cDBGSerialization.enum_from_str(cdbg_format),
                                                  self.compactor.cdbg.__smartptr__())
                self.processor.register_listener(writer)
                self.writers.append(writer)
                           
    def execute(self, args):
        for sample, prefix in iter_fastx_inputs(args.inputs, args.pairing_mode):
            self.processor.process(*sample)

    def teardown(self):
        pass

