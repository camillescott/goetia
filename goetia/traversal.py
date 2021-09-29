#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : traversal.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.10.2019

from goetia import libgoetia
from goetia.utils import copy_attrs
from cppyy.gbl import std
from cppyy import gbl

UnitigWalker = libgoetia.UnitigWalker
STATES = libgoetia.TraversalState
