#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : runner.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 11.03.2020

import sys

import blessings


class CommandRunner:

    def __init__(self, parser, *args, description='', **kwargs):
        self.parser = parser
        self.description = description
        self.term = blessings.Terminal()

    def postprocess_args(self, args):
        pass

    def setup(self, args):
        pass

    def execute(self, args):
        raise NotImplementedError

    def teardown(self):
        pass

    def run(self, args):
        try:
            desc = self.description.format(term=self.term)
        except:
            desc = ''
        finally:
            print(desc, file=sys.stderr)

        self.postprocess_args(args)
        self.setup(args)
        self.execute(args)
        self.teardown()
