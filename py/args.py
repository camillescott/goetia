# This file is adapted from khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org

import sys
import argparse
import math
import textwrap
from argparse import _VersionAction
from collections import namedtuple

from boink import __version__
from boink.serialization import cDBGSerialization
from boink.metadata import CUR_TIME
from boink.parsing import PAIRING_MODES
from boink import types
from boink import libboink


DEFAULT_K = 31
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = 1e8
DEFAULT_N_THREADS = 1


def memory_setting(label):
    """
    Parse user-supplied memory setting.

    Converts strings into floats, representing a number of bytes. Supports the
    following notations.
      - raw integers: 1, 1000, 1000000000
      - scientific notation: 1, 1e3, 1e9
      - "common" notation: 1, 1K, 1G

    Suffixes supported: K/k, M/m, G/g, T/t. Do not include a trailing B/b.
    """
    suffixes = {
        'K': 1000.0,
        'M': 1000.0 ** 2,
        'G': 1000.0 ** 3,
        'T': 1000.0 ** 4,
    }
    try:
        mem = float(label)
        return mem
    except ValueError:
        prefix = label[:-1]
        suffix = label[-1:].upper()
        if suffix not in suffixes.keys():
            raise ValueError('cannot parse memory setting "{}"'.format(label))
        try:
            multiplier = float(prefix)
            return multiplier * suffixes[suffix]
        except ValueError:
            raise ValueError('cannot parse memory setting "{}"'.format(label))


class _VersionStdErrAction(_VersionAction):
    # pylint: disable=too-few-public-methods, protected-access
    """Force output to StdErr."""

    def __call__(self, parser, namespace, values, option_string=None):
        # have to call info() directly as the version action exits
        # which means parse_args() does not get a chance to run
        version = self.version
        if version is None:
            version = parser.version
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), sys.stderr)
        parser.exit()


class ComboFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    """Both ArgumentDefaults and RawDescription formatters."""

    pass


class BoinkArgumentParser(argparse.ArgumentParser):
    """Specialize ArgumentParser with boink defaults.

    Take care of common arguments and setup printing of citation information.
    """

    def __init__(self, formatter_class=ComboFormatter,
                 **kwargs):
        super(BoinkArgumentParser, self).__init__(
            formatter_class=formatter_class, add_help=False, **kwargs)

        self.add_argument('--version', action=_VersionStdErrAction,
                          version='boink {v}'.format(v=__version__))
        self.add_argument('-h', '--help',
                          default=argparse.SUPPRESS,
                          help='show this help message and exit')


def print_boink_intro():
    print('*' * 20, '*', sep='\n', file=sys.stderr)
    print('*    BOINK v{0}'.format(__version__), file=sys.stderr)
    print('*' * 20, '*', sep='\n', file=sys.stderr)


def build_dBG_args(descr=None, epilog=None, parser=None):
    """Build an ArgumentParser with args for bloom filter based scripts."""
    expert_help = '--help-expert' in sys.argv
    if expert_help:
        sys.argv.append('--help')

    if parser is None:
        parser = BoinkArgumentParser(description=descr, epilog=epilog)

    DBG_TYPES    = [t[0].__name__ for t in types.storage_types]
    DEFAULT_ARGS = {t[0].__name__ : t[1] for t in types.storage_types}

    class dBGArgAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            pass

    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K,
                        help='k-mer size to use')

    help = ('number of tables to use in k-mer countgraph' if expert_help
            else argparse.SUPPRESS)
    parser.add_argument('--storage-type',
                        default=DBG_TYPES[0],
                        choices=DBG_TYPES)
    parser.add_argument('--n_tables', '-N', type=int,
                        default=DEFAULT_N_TABLES,
                        help=help)

    parser.add_argument('-U', '--unique-kmers', type=float, default=0,
                        help='approximate number of unique kmers in the input'
                             ' set')
    parser.add_argument('--fp-rate', type=float, default=None,
                        help="Override the automatic FP rate setting for the"
                        " current script")

    group = parser.add_mutually_exclusive_group()
    help = ('upper bound on tablesize to use; overrides --max-memory-usage/-M'
            if expert_help else argparse.SUPPRESS)
    group.add_argument('--max-tablesize', '-x', type=float,
                       default=DEFAULT_MAX_TABLESIZE,
                       help=help)
    group.add_argument('-M', '--max-memory-usage', type=memory_setting,
                       help='maximum amount of memory to use for data ' +
                       'structure')

    group.add_argument('--save-dbg', metavar="filename", default=None,
                        help='save the k-mer countgraph to disk after all '
                        'reads are loaded.')

    return parser


def print_dBG_args(args):
    print('* dBG will be order', args.ksize, file=sys.stderr)
    print('* dBG will have underlying storage of', args.storage_type, file=sys.stderr)
    print('*','*' * 10, '*', sep='\n', file=sys.stderr)


def add_output_interval_args(parser):
    parser.add_argument('--fine-interval', type=int, default=libboink.DEFAULT_INTERVALS.FINE)
    parser.add_argument('--medium-interval', type=int, default=libboink.DEFAULT_INTERVALS.MEDIUM)
    parser.add_argument('--coarse-interval', type=int, default=libboink.DEFAULT_INTERVALS.COARSE)


def print_interval_settings(args):
    print('* FINE output interval:', args.fine_interval, file=sys.stderr)
    print('* MEDIUM output interval:', args.medium_interval, file=sys.stderr)
    print('* COARSE output interval:', args.coarse_interval, file=sys.stderr)
    print('*', '*' * 10, '*', sep='\n', file=sys.stderr)


def add_save_cDBG_args(parser):
    default_prefix = 'boink.build-cdbg.' + CUR_TIME
    parser.default_prefix = default_prefix

    parser.add_argument('--results-dir',
                        default=default_prefix)

    parser.add_argument('--normalize',
                        type=int,
                        nargs='?',
                        const=10)

    parser.add_argument('--save-cdbg',
                        metavar='PREFIX.<format>',
                        nargs='?',
                        const='boink.cdbg.graph')
    parser.add_argument('--save-cdbg-format',
                        nargs='+',
                        choices=cDBGSerialization.FORMATS,
                        default=['gfa1'])

    parser.add_argument('--track-cdbg-stats',
                        metavar='FILE_NAME.csv',
                        nargs='?',
                        const='boink.cdbg.stats.csv')

    parser.add_argument('--track-cdbg-history',
                        metavar='FILENAME.graphml',
                        nargs='?',
                        const='boink.cdbg.history.graphml')

    parser.add_argument('--track-cdbg-components',
                        metavar='FILE_NAME.csv',
                        nargs='?',
                        const='boink.cdbg.components.csv')
    parser.add_argument('--component-sample-size',
                        type=int,
                        default=10000)

    parser.add_argument('--track-cdbg-unitig-bp',
                        metavar='FILENAME.csv',
                        nargs='?',
                        const='boink.cdbg.unitigs.bp.csv')

    parser.add_argument('--unitig-bp-bins',
                        nargs='+',
                        type=int)

    parser.add_argument('--validate',
                        metavar='FILENAME.csv',
                        nargs='?',
                        const='boink.cdbg.validation.csv')

    return parser


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


def add_prometheus_args(parser):
    parser.add_argument('--port', default=None,
                        help='Port to expose prometheus metrics.')
    return parser


def print_prometheus_args(args):
    print('* Exposing prometheus metrics on port', args.port, file=sys.stderr)
    print('*', '*' * 10, '*', sep='\n', file=sys.stderr)


def add_pairing_args(parser):
    """Common pairing mode argument."""
    parser.add_argument('--pairing-mode', default='interleaved',
                        choices=PAIRING_MODES,
                        help='How to interpret read pairing. With `single`, '\
                             'reads will be parsed as singletons, regardless'\
                             ' of pairing or file order. With `interleaved`,'\
                             ' each file will be assumed to be interleaved '\
                             'and paired, with singletons allowed to be mixed'\
                             ' in. With `split`, it will be assumed that each'\
                             ' group of two files in the input list are '\
                             'as (LEFT, RIGHT), ...')
    return parser
