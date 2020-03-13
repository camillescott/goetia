#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : cli.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

import argparse
import sys

from goetia import __version__
from goetia.processors import DEFAULT_INTERVALS


class _VersionStdErrAction(argparse._VersionAction):
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


class GoetiaArgumentParser(argparse.ArgumentParser):
    """Specialize ArgumentParser with goetia defaults.
    """

    def __init__(self, formatter_class=ComboFormatter,
                 **kwargs):
        super(GoetiaArgumentParser, self).__init__(
            formatter_class=formatter_class, **kwargs)

        self.add_argument('--version', action=_VersionStdErrAction,
                          version='goetia {v}'.format(v=__version__))
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help(sys.stderr)
        sys.exit(2)


def print_goetia_intro():
    print('*' * 20, '*', sep='\n', file=sys.stderr)
    print('*    GOETIA v{0}'.format(__version__), file=sys.stderr)
    print('*' * 20, '*', sep='\n', file=sys.stderr)


def get_output_interval_args(parser):
    group = parser.add_argument_group('reporting')
    group.add_argument('--fine-interval', type=int, default=DEFAULT_INTERVALS.FINE)
    group.add_argument('--medium-interval', type=int, default=DEFAULT_INTERVALS.MEDIUM)
    group.add_argument('--coarse-interval', type=int, default=DEFAULT_INTERVALS.COARSE)
    return group


def print_interval_settings(args):
    print('* FINE output interval:', args.fine_interval, file=sys.stderr)
    print('* MEDIUM output interval:', args.medium_interval, file=sys.stderr)
    print('* COARSE output interval:', args.coarse_interval, file=sys.stderr)
    print('*', '*' * 10, '*', sep='\n', file=sys.stderr)

