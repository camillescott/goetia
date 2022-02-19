#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : cdbg.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 14.10.2019

import json
import sys

from boltons.iterutils import pairwise
import curio

from goetia import libgoetia
from goetia.serialization import cDBGSerialization
from goetia.timer import time_block

cDBG = libgoetia.cDBG
StreamingCompactor = libgoetia.StreamingCompactor


def quick_compactor(K):
    from goetia.dbg import dBG
    from goetia.storage import PHMapStorage
    from goetia.hashing import FwdLemireShifter

    storage = PHMapStorage.build()
    hasher = FwdLemireShifter(K)
    dbg = dBG[type(storage), type(hasher)].build(storage, hasher)
    compactor_t = StreamingCompactor[type(dbg)]
    compactor = compactor_t.Compactor.build(dbg)

    return compactor, compactor_t


async def compute_connected_component_callback(msg,
                                               cdbg_type,
                                               cdbg,
                                               out_stream,
                                               sample_size,
                                               verbose=None):

    if verbose:
        verbose.message(f'Computing sample of connected component sizes...')

    with time_block() as elapsed:
        result = cdbg_type.compute_connected_component_metrics(cdbg,
                                                               sample_size)
        n_comps, min_comp, max_comp, sizes = result

        data = {'t': msg.t,
                'sample_name': msg.sample_name,
                'n_components': n_comps,
                'max': max_comp,
                'min': min_comp,
                'sizes': list(sizes)}

    if verbose:
        verbose.message(f'Found {n_comps} components, minimum size {min_comp}, '\
                        f'maximum size {max_comp} (took {elapsed}).')

    await out_stream.write(data)


async def compute_unitig_fragmentation_callback(msg,
                                                cdbg_type,
                                                cdbg,
                                                out_stream,
                                                bins,
                                                verbose=None):

    if verbose:
        verbose.message(f'Computing unitig length hist (bins={bins})...')

    counts = cdbg_type.compute_unitig_fragmentation(cdbg,
                                                    bins)
    data = {'t': msg.t,
            'sample_name': msg.sample_name}
    for i in range(len(bins) - 1):
        bin_start, bin_end = bins[i], bins[i+1]
        data[f'[{bin_start},{bin_end})'] = counts[i]
    data[f'[{bin_end},Inf)'] = counts[-1]

    if verbose:
        verbose.message(f'Bin counts: {data}')

    await out_stream.write(data)


async def write_cdbg_metrics_callback(msg,
                                      compactor,
                                      out_stream,
                                      verbose=None):

    report = compactor.get_report()
    data = {'t': msg.t,
            'seq_t': msg.sequence,
            'rt_elapsed_interval': msg.seconds_elapsed_interval,
            'rt_elapsed_total': msg.seconds_elapsed_total,
            'sample_name': msg.sample_name,
            'n_full': report.n_full,
            'n_tips': report.n_tips,
            'n_islands': report.n_islands,
            'n_trivial': report.n_trivial,
            'n_circular': report.n_circular,
            'n_loops': report.n_loops,
            'n_dnodes': report.n_dnodes,
            'n_unodes': report.n_unodes,
            'n_tags': report.n_tags,
            'n_updates': report.n_updates,
            'n_splits': report.n_splits,
            'n_merges': report.n_merges,
            'n_extends': report.n_extends,
            'n_clips': report.n_clips,
            'n_deletes': report.n_deletes,
            'n_circular_merges': report.n_circular_merges,
            'n_unique_kmers': report.n_unique,
            'estimated_fp': report.estimated_fp}

    await out_stream.write(data)


async def write_cdbg_callback(msg,
                              cdbg,
                              out_file_prefix,
                              out_format,
                              time_component=None,
                              verbose=None):

    if time_component is None:
        out_file_name = f'{out_file_prefix}.{out_format}'
    else:
        out_file_name = f'{out_file_prefix}.{msg.t}.{out_format}'

    if verbose:
        verbose.message(f'Writing cDBG to {out_file_name}...')

    cdbg.write(out_file_name,
                cDBGSerialization.enum_from_str(out_format))
    
    if verbose:
        verbose.message(f'Finished writing cDBG.', file=sys.stderr)


async def validate_cdbg_callback(msg,
                                 cdbg,
                                 out_file_name):
    cdbg.validate(out_file_name)
