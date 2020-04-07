#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : saturation.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 03.03.2020

from abc import ABC, abstractmethod
from collections import deque
from statistics import median

import numpy as np


class SaturationTracker(ABC):

    def __init__(self, norm_tick_length, window_size, dfunc):
        self.window_size      = window_size
        self.norm_tick_length = norm_tick_length
        self.saturated        = False
        self.distances        = [0.0]
        self.read_times       = [0]
        self.norm_times       = [0]
        self.dfunc            = dfunc
        self.prev_sig         = None

    @abstractmethod
    def push(self, new_sig, new_read_time, new_norm_time):
        pass

    def find_norm_index(self):
        if len(self.norm_times) < 2:
            print('len norm times < 2')
            return None, None

        i = 2
        norm_sum = 0
        try:
            while True:
                norm_sum = norm_sum + (self.norm_times[-i + 1] - self.norm_times[-i])
                if norm_sum > self.norm_tick_length:
                    return -i, norm_sum
                i += 1
        except IndexError as e:
            print(e, i)
            return None, None


class MedianDistanceSaturation(SaturationTracker):

    def __init__(self, norm_tick_length, window_size, distance_cutoff, dfunc):
        self.distance_cutoff = distance_cutoff
        self.medians         = []

        super().__init__(norm_tick_length, window_size, dfunc)

    def push(self, new_sig, new_read_time, new_norm_time):
        if self.prev_sig:
            self.distances.append(self.dfunc(self.prev_sig, new_sig))
            self.read_times.append(new_read_time)
            self.norm_times.append(new_norm_time)
            retval = self.distances[-1], self.read_times[-1], np.NaN
        else:
            retval = None, None, None
        self.prev_sig = new_sig

        window_index, norm_sum = self.find_norm_index()
        if window_index is not None:
            med = median(self.distances[window_index:])
            self.medians.append(med)
            retval = self.distances[-1], self.read_times[-1], med
            if len(self.medians) >= self.window_size and \
                    all((d > self.distance_cutoff for d in self.medians[-self.window_size:])):
                self.saturated = True

        return retval


SaturationPolicies = {'median.distance': MedianDistanceSaturation}
