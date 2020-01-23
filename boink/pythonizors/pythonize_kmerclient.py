#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : pythonize_kmerclient.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 22.01.2020


from cppyy import gbl


def pythonize_boink_kmers(klass, name):

    if name == 'KmerClient':

        klass.K = property(klass.K)
