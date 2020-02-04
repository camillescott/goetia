#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : traversal.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.10.2019

from boink import libboink
from boink.utils import copy_attrs
from cppyy.gbl import std
from cppyy import gbl

dBGWalker = libboink.dBGWalker
STATES = libboink.TraversalState
