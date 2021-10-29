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

from boltons.iterutils import windowed_iter
import ijson
import numpy as np
from scipy.spatial.distance import cosine
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

    def similarity(self, other, metric=cosine):
        return metric(self.sketch, other.sketch)
    
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
                'name': self.name(),
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


def distance(sig_a, sig_b, metric=cosine):
    if isinstance(sig_a, SourmashSignature):
        return sig_a.similarity(sig_b)
    elif isinstance(sig_a, DraffSignature):
        return sig_a.similarity(sig_b, metric=metric)
    else:
        raise TypeError(f'Not a support signature type: {type(sig_a)}.')


def find_rolling_distances(sigs, dmetrics=[cosine], window_sizes = [2,4,6,8,10]):
    '''

    '''
    max_window = max(window_sizes)
    window_sizes.sort()
    times, distances, freqs, metrics = [], [], [], []
    window_freqs = {}

    for i, window in enumerate(windowed_iter(sigs, max_window)):
        for sub_window_size in window_sizes:
            sig_a, sig_b = window[0], window[sub_window_size - 1]
            t_a = int(sig_a.name().split(':')[1])
            t_b = int(sig_b.name().split(':')[1])
            if i == 0:
                freq = round(t_b - t_a, -4)
                window_freqs[sub_window_size] = freq
            else:
                freq = window_freqs[sub_window_size]

            for metric in dmetrics:
                times.append(int(t_b))
                distances.append(distance(sig_a, sig_b, metric=metric))
                metrics.append(metric.__name__)
                freqs.append(int(freq))

            #print(window[0], window[sub_window_size - 1])

    for sub_window_size in window_sizes[:-1]:
        for sub_window in windowed_iter(window[1:], sub_window_size):
            #print(sub_window[0], sub_window[-1])
            sig_a, sig_b = sub_window[0], sub_window[-1]
            t_a = int(sig_a.name().split(':')[1])
            t_b = int(sig_b.name().split(':')[1])
            
            for metric in dmetrics:
                times.append(int(t_b))
                distances.append(distance(sig_a, sig_b, metric=metric))
                metrics.append(metric.__name__)
                freqs.append(window_freqs[sub_window_size])

    df = pd.DataFrame({'time': times,
                       'distance': distances,
                       'freq': freqs,
                       'metric': metrics})
    return df


def find_distances_from_ref(sigs, ref_sig, dmetrics=[cosine], cutoff=1.0):
    times, distances, metrics = [], [], []
    for sig in sigs:
        sample_t = int(sig.name().split(':')[1])

        for metric in dmetrics:
            times.append(sample_t)
            distances.append(distance(sig, ref_sig, metric=metric))
            metrics.append(metric.__name__)
        
    df = pd.DataFrame({'time': times,
                       'distance': distances,
                       'metric': metrics})
    return df
