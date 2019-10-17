#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : cli.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

import argparse
import sys

from boink import __version__, libboink


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


class BoinkArgumentParser(argparse.ArgumentParser):
    """Specialize ArgumentParser with boink defaults.
    """

    def __init__(self, formatter_class=ComboFormatter,
                 **kwargs):
        super(BoinkArgumentParser, self).__init__(
            formatter_class=formatter_class, **kwargs)

        self.add_argument('--version', action=_VersionStdErrAction,
                          version='boink {v}'.format(v=__version__))
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help(sys.stderr)
        sys.exit(2)



def print_boink_intro():
    print('*' * 20, '*', sep='\n', file=sys.stderr)
    print('*    BOINK v{0}'.format(__version__), file=sys.stderr)
    print('*' * 20, '*', sep='\n', file=sys.stderr)



def get_output_interval_args(parser):
    group = parser.add_argument_group('reporting')
    group.add_argument('--fine-interval', type=int, default=libboink.DEFAULT_INTERVALS.FINE)
    group.add_argument('--medium-interval', type=int, default=libboink.DEFAULT_INTERVALS.MEDIUM)
    group.add_argument('--coarse-interval', type=int, default=libboink.DEFAULT_INTERVALS.COARSE)
    return group


def print_interval_settings(args):
    print('* FINE output interval:', args.fine_interval, file=sys.stderr)
    print('* MEDIUM output interval:', args.medium_interval, file=sys.stderr)
    print('* COARSE output interval:', args.coarse_interval, file=sys.stderr)
    print('*', '*' * 10, '*', sep='\n', file=sys.stderr)


class BoinkRunner:

    def __init__(self):
        self.parser = BoinkArgumentParser()
        self.parser.set_defaults(func = lambda _: self.parser.print_help())
        self.commands = self.parser.add_subparsers()

    def add_command(self, name, runner_type):
        cmd_parser = self.commands.add_parser(name)
        cmd_runner = runner_type(cmd_parser)
        cmd_parser.set_defaults(func=cmd_runner.run)
        return cmd_parser

    def run(self):
        args = self.parser.parse_args()
        return args.func(args)


class CommandRunner:

    def __init__(self, parser, *args, **kwargs):
        self.parser = parser

    def postprocess_args(self, args):
        pass

    def setup(self, args):
        pass

    def execute(self, args):
        raise NotImplementedError

    def teardown(self):
        pass

    def run(self, args):
        self.postprocess_args(args)
        self.setup(args)
        self.execute(args)
        self.teardown()
