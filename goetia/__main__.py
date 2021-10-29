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
from goetia.cli.diginorm import DiginormFilterRunner
from goetia.cli.solid_filter import SolidFilterRunner
from goetia.cli.sourmash_stream import SourmashRunner
from goetia.cli.draff_stream import DraffRunner
from goetia.cli.merge_paired import PairedMergeRunner


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


    # parent for `goetia cdbg`
    cdbg = commands.add_parser('cdbg')
    cdbg_commands = cdbg.add_subparsers()
    
    # `goetia cdbg`
    cdbg_build_parser = cdbg_commands.add_parser('build')
    cdbg_build_command = cDBGRunner(cdbg_build_parser)
    cdbg_build_parser.set_defaults(func=about(cdbg_build_command.run))

    # parent for `goetia filter`
    filters = commands.add_parser('filter')
    filter_commands = filters.add_subparsers()

    # `goetia filter solid`
    solid_filter_parser = filter_commands.add_parser('solid')
    solid_filter_command = SolidFilterRunner(solid_filter_parser)
    solid_filter_parser.set_defaults(func=about(solid_filter_command.run))

    # `goetia filter diginorm`
    diginorm_filter_parser = filter_commands.add_parser('diginorm')
    diginorm_filter_command = DiginormFilterRunner(diginorm_filter_parser)
    diginorm_filter_parser.set_defaults(func=about(diginorm_filter_command.run))

    # parent for `goetia utils`
    utils = commands.add_parser('utils')
    utils_commands = utils.add_subparsers()

    # `goetia utils merge-paired`
    merge_split_parser = utils_commands.add_parser('merge-paired')
    merge_split_command = PairedMergeRunner(merge_split_parser)
    merge_split_parser.set_defaults(func=about(merge_split_command.run))

    # run the command
    args = parser.parse_args()
    return args.func(args)

