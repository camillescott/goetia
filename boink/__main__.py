#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : __main__.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

from boink.cli import BoinkRunner
from boink.cdbg import cDBGRunner


def main():
    runner = BoinkRunner()
    runner.add_command('cdbg', cDBGRunner)

    return runner.run()
