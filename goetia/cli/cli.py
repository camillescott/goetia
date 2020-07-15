#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : cli.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

import argparse
import os
import sys

from goetia import __version__


def print_goetia_intro():
    print('*' * 20, '*', sep='\n', file=sys.stderr)
    print('*    GOETIA v{0}'.format(__version__), file=sys.stderr)
    print('*' * 20, '*', sep='\n', file=sys.stderr)


def format_filenames(file_names):
    return '+'.join([os.path.basename(name) for name in file_names])


class TextBlock:

    def __init__(self, text):
        self.text = text.split('\n')
        self.height = len(self.text)
        self.width = max((len(line) for line in self.text))

    def __iter__(self):
        yield from self.text

    def draw(self, term, x, shift_y=None):
        out = term.stream
        
        if shift_y:
            buf = term.move_down * shift_y
        else:
            buf = ''

        for line in self.text:
            buf += term.move_x(x)
            buf += line
            buf += term.move_down
        
        out.write(buf)



