#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : saturation.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 03.03.2020

from collections import deque
from statistics import median, mean

import numpy as np


class SlidingWindow:

    def __init__(self, window_size, func, uses_time=False):
        if window_size < 2:
            raise TypeError('window_size must be at least 2')
        self.window_size = window_size
        self.func = func
        self.values = []
        self.uses_time = uses_time
    
    def reset(self):
        self.values = []
    
    def push(self, value):
        self.values.append(self._unpack(value))
        
        if len(self.values) >= self.window_size:
            stat = self.func([v for v, t in self.values[-self.window_size:]]) if not self.uses_time else \
                   self.func(self.values[-self.window_size:])
        else:
            stat = np.NaN
        
        _, back_time = self.values[-1]
        if np.isnan(back_time):
            back_time = len(self.values) - 1

        return stat, back_time
    
    @staticmethod
    def _unpack(item):
        try:
            value, time = item
        except TypeError:
            return item, np.NaN
        else:
            return value, time


class SlidingCutoff:

    def __init__(self, window_size, window_func, cutoff_func):
        self.value_window = SlidingWindow(window_size, window_func)
        self.window_value = np.NaN
        self.cutoff_window = SlidingCutoff(window_size, cutoff_func)
        self.cutoff_reached = False
        self.time = 0

    def push(self, value):
        self.window_value, self.time = self.value_window.push(value)
        if not np.isnan(self.window_value):
            self.cutoff_reached, _ = self.cutoff_window.push((self.window_value, self.time))
            return self.cutoff_reached, self.window_value, self.time
        else:
            return False, self.window_value, self.time


def all_cutoff(cutoff):
    def func(values):
        return all((v > cutoff for v in values))
    return func


def median_cutoff(cutoff):
    def func(values):
        return median(values) > cutoff


def normalized_mean(values):
    m = mean((v for v, t in values))
    n = values[-1][-1] - values[0][-1]

    return m / n


cutoff_functions = { 'all': all_cutoff,
                     'median': median_cutoff }


smoothing_functions = { 'mean': np.mean,
                        'median': median,
                        'stddev': np.std }