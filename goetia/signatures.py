#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 15.10.2019

import hashlib
import json
import typing

from sourmash import SourmashSignature

from goetia import __version__

class DraffSignature:

    def __init__(self, sketch, name='', filename='', license='CC0'):
        self._sketch = sketch.to_numpy()

        if name:
            self._name = name
        if filename:
            self._filename = filename

        self.license = license
        self.W = sketch.W
        self.K = sketch.K
    
    def md5sum(self):
        return hashlib.md5(self._sketch.data).hexdigest()
    
    @property
    def name(self):
        if self._name:
            return self._name
        if self._filename:
            return self._filename
        else:
            return self.md5sum()[:8]

    @property
    def sketch(self):
        return self._sketch
    
    @property
    def size(self):
        return len(self._sketch)
    
    def to_dict(self):
        return {'sketch': self._sketch.to_list(),
                'W': self.W,
                'K': self.K,
                'size': self.size,
                'version': __version__,
                'license': self.license,
                'md5sum': self.md5sum()}

    def save(self, stream: typing.TextIO) -> None:
        """Save the signature to disk.

        Args:
            stream (file): File handle to save to.
            name (str): Name of the signature.
        """
        data = [self.to_dict()]

        json.dump(data, stream)
