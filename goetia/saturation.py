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

    def __init__(self, window_size, func, uses_time=False, history=0):
        if window_size < 2:
            raise TypeError('window_size must be at least 2')
        self.window_size = window_size
        self.func = func
        self.items = []
        self.uses_time = uses_time

        if history > 0 and history < window_size:
            history = window_size
        self.history = history
    
    @property
    def name(self):
        return f'SlidingWindow({self.func.__name__}, L={self.window_size})'
    
    def reset(self):
        self.items = []
    
    def times(self):
        for _, time in self.items:
            yield time
    
    def values(self):
        for value, _ in self.items:
            yield value
    
    def push(self, value):
        self.items.append(self._unpack(value))
        
        if len(self.items) >= self.window_size:
            stat = self.func([v for v, t in self.items[-self.window_size:]]) if not self.uses_time else \
                   self.func(self.items[-self.window_size:])
        else:
            stat = np.NaN
        
        _, back_time = self.items[-1]
        if np.isnan(back_time):
            back_time = len(self.items) - 1

        if self.history > 0:
            self.items = self.items[-self.history:]

        return stat, back_time
    
    @staticmethod
    def _unpack(item):
        try:
            value, time = item
        except TypeError:
            return item, np.NaN
        else:
            return value, time


class RollingPairwise(SlidingWindow):

    def __init__(self, dfunc, **kwargs):
        super().__init__(2, dfunc, **kwargs)


class SlidingCutoff:

    def __init__(self, window_size, window_func, cutoff_func):
        if isinstance(window_func, SlidingWindow):
            self.value_window = window_func
        else:
            self.value_window = SlidingWindow(window_size, window_func)
        self.window_value = np.NaN

        if isinstance(cutoff_func, SlidingWindow):
            self.cutoff_window = cutoff_func
        else:
            self.cutoff_window = SlidingWindow(window_size, cutoff_func)
        self.cutoff_reached = False

        self.time = 0
    
    @property
    def name(self):
        return f'SlidingCutoff(smoothing={self.value_window.name}, cutoff={self.cutoff_window.name})'
    
    @property
    def cutoff_reached(self):
        return self._cutoff_reached
    
    @cutoff_reached.setter
    def cutoff_reached(self, reached):
        if np.isnan(reached):
            self._cutoff_reached = False
        else:
            self._cutoff_reached = bool(reached)

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
    return func


def normalized_mean(items):
    m = mean((v for v, t in items))
    n = items[-1][-1] - items[0][-1]

    return m / n


cutoff_functions = { 'all': all_cutoff,
                     'median': median_cutoff }


smoothing_functions = { 'mean': np.mean,
                        'median': median,
                        'stddev': np.std }
