#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : __main__.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

from goetia.cli import GoetiaRunner
from goetia.cdbg import cDBGRunner
from goetia.signatures import SourmashRunner


def main():
    runner = GoetiaRunner()
    runner.add_command('cdbg', cDBGRunner)
    runner.add_command('sourmash', SourmashRunner)

    return runner.run()
