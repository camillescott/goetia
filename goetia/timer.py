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


@contextmanager
def measure_time(fp=sys.stdout, *args, **kwargs):
    start = time.perf_counter()
    yield
    end = time.perf_counter()
    elapsed = end - start
    print(f'Elapsed: {elapsed:0.4f}', file=fp)
