#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : runner.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 11.03.2020

from goetia.cli.args import GoetiaArgumentParser


class GoetiaRunner:

    def __init__(self):
        self.parser = GoetiaArgumentParser()
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
