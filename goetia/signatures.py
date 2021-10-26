#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : signatures.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 15.10.2019

import decimal
import hashlib
import json
import typing

import ijson
import numpy as np
from sourmash import SourmashSignature
from sourmash._lowlevel import ffi, lib
from sourmash.utils import RustObject, rustcall, decode_str

from goetia import __version__


class DraffSignature:

    def __init__(self, sketch, name='', filename='', license='CC0',
                 W=None, K=None, version=__version__):
        if not isinstance(sketch, np.ndarray):
            self._sketch = sketch.to_numpy()
        else:
            self._sketch = sketch

        if name:
            self._name = name
        if filename:
            self._filename = filename

        self.version = version
        self.license = license
        self.W = sketch.W if W is None else W
        self.K = sketch.K if K is None else K
    
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
        return {'sketch': self._sketch.tolist(),
                'W': self.W,
                'K': self.K,
                'size': self.size,
                'name': self.name,
                'version': self.version,
                'license': self.license,
                'md5sum': self.md5sum()}

    @classmethod
    def from_dict(cls, data):
        sig = cls(np.array(data['sketch']),
                  name = data['name'],
                  license = data['license'],
                  version = data['version'],
                  W = int(data['W']),
                  K = int(data['K']))
        return sig

    def save(self, stream: typing.TextIO) -> None:
        """Save the signature to disk.

        Args:
            stream (file): File handle to save to.
            name (str): Name of the signature.
        """
        data = [self.to_dict()]

        json.dump(data, stream)


class DecimalEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, decimal.Decimal):
            return float(o)
        return super(DecimalEncoder, self).default(o)


def load_draff_stream(fp):
    '''Iteratively parser a JSON file of draff
    signatures from the given file.

    Args:
        fp (file): File handle to parse from.
    '''

    backend = ijson.get_backend('yajl2_c')
    for signature in backend.items(fp, 'item'):
        yield DraffSignature.from_dict(signature)


def load_sourmash_stream(fp):
    '''Iteratively parse a JSON file of sourmash
    signatures from the given file.

    Args:
        fp (file): File handle to parse from.
    '''

    backend = ijson.get_backend('yajl2_c')
    for signature in backend.items(fp, 'item'):
        data = json.dumps([signature], cls=DecimalEncoder).encode('utf-8')

        size = ffi.new("uintptr_t *")
        ptr = rustcall(lib.signatures_load_buffer, data, len(data), False, 0, ffi.NULL, size)
        size = ffi.unpack(size, 1)[0]
        sigs = []
        for i in range(size):
            sigs.append(SourmashSignature._from_objptr(ptr[i]))
        yield sigs[0]
