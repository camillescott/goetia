#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : __main__.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

import sys

from goetia import splash
from goetia.cli.args import GoetiaArgumentParser
from goetia.cli.cdbg_stream import cDBGRunner
from goetia.cli.solid_filter import SolidFilterRunner
from goetia.cli.sourmash_stream import SourmashRunner
from goetia.cli.draff_stream import DraffRunner


def about(func):
    def command(*args, **kwargs):
        print(splash, file=sys.stderr)
        return func(*args, **kwargs)
    return command


def main():
    parser = GoetiaArgumentParser(
        description=splash
    )
    parser.set_defaults(func = lambda _: self.parser.print_help())
    commands = parser.add_subparsers()
    
    # parent for `goetia sketch`
    sketch = commands.add_parser('sketch')
    sketch_commands = sketch.add_subparsers()

    # `goetia sketch sourmash`
    sketch_sourmash_parser = sketch_commands.add_parser('sourmash')
    sketch_sourmash_command = SourmashRunner(sketch_sourmash_parser)
    sketch_sourmash_parser.set_defaults(func=about(sketch_sourmash_command.run))

    # `goetia sketch draff`
    sketch_draff_parser = sketch_commands.add_parser('draff')
    sketch_draff_command = DraffRunner(sketch_draff_parser)
    sketch_draff_parser.set_defaults(func=about(sketch_draff_command.run))

    # `goetia cdbg`
    cdbg_parser = commands.add_parser('cdbg')
    cdbg_command = cDBGRunner(cdbg_parser)
    cdbg_parser.set_defaults(func=about(cdbg_command.run))

    # `goetia solid-filter`
    solid_filter_parser = commands.add_parser('solid-filter')
    solid_filter_command = SolidFilterRunner(solid_filter_parser)
    solid_filter_parser.set_defaults(func=about(solid_filter_command.run))

    # run the command
    args = parser.parse_args()
    return args.func(args)

