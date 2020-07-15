#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : __main__.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

from goetia.cli.runner import GoetiaRunner
from goetia.cli.cdbg_stream import cDBGRunner
from goetia.cli.solid_filter import SolidFilterRunner
from goetia.cli.sourmash_stream import SourmashRunner


def main():
    runner = GoetiaRunner()
    runner.add_command('cdbg', cDBGRunner)
    runner.add_command('sourmash', SourmashRunner)
    runner.add_command('solid-filter', SolidFilterRunner)

    return runner.run()
