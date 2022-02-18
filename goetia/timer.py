#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : timer.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 12.03.2021

from contextlib import contextmanager
import sys
import time


class Elapsed:
    def __init__(self):
        self.start = time.perf_counter()
        self.end = None
        self.elapsed = -1

    def stop(self):
        self.end = time.perf_counter()
        self.elapsed = self.end - self.start

    def __str__(self):
        return f'{self.elapsed:0.2f}s'

    def __float__(self):
        return self.elapsed


@contextmanager
def time_block():
    e = Elapsed()
    yield e
    e.stop()
