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

cDBG = libgoetia.cDBG
StreamingCompactor = libgoetia.StreamingCompactor


async def compute_connected_component_callback(msg,
                                               cdbg_type,
                                               cdbg,
                                               out_file,
                                               sample_size,
                                               verbose=False):

    if verbose:
        print(f'Computing sample of connected component sizes...',
              file=sys.stderr)

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
        print(f'Found {n_comps} components, of minimum size {min_comp}, '\
              f'maximum size {max_comp}.',
              file=sys.stderr)

    async with curio.aopen(out_file, 'a') as fp:
        if await fp.tell() != 0:
            await fp.write(',\n')
        else:
            await fp.write('[\n')
        await fp.write(json.dumps(data))


async def compute_unitig_fragmentation_callback(msg,
                                                cdbg_type,
                                                cdbg,
                                                out_file,
                                                bins,
                                                verbose=False):

    if verbose:
        print(f'Computing unitig length hist (bins={bins})...',
              file=sys.stderr)

    counts = cdbg_type.compute_unitig_fragmentation(cdbg,
                                                    bins)
    data = {'t': msg.t,
            'sample_name': msg.sample_name}
    for i in range(len(bins) - 1):
        bin_start, bin_end = bins[i], bins[i+1]
        data[f'[{bin_start},{bin_end})'] = counts[i]
    data[f'[{bin_end},Inf)'] = counts[-1]

    if verbose:
        print(f'Bin counts: {data}',
              file=sys.stderr)

    async with curio.aopen(out_file, 'a') as fp:
        if await fp.tell() != 0:
            await fp.write(',\n')
        else:
            await fp.write('[\n')
        await fp.write(json.dumps(data))


async def write_cdbg_metrics_callback(msg,
                                      compactor,
                                      out_file,
                                      verbose=False):

    report = compactor.get_report()
    data = {'t': msg.t,
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
    
    if verbose:
        print(f'Writing metrics: {report.n_dnodes} decision nodes, {report.n_updates} updates...',
              file=sys.stderr)

    async with curio.aopen(out_file, 'a') as fp:
        if await fp.tell() != 0:
            await fp.write(',\n')
        else:
            await fp.write('[\n')
        await fp.write(json.dumps(data))


async def write_cdbg_callback(msg,
                              cdbg,
                              out_file_prefix,
                              out_format,
                              verbose=False):

    out_file_name = f'{out_file_prefix}.{msg.t}.{out_format}'

    if verbose:
        print(f'Writing cDBG to {out_file_name}...', file=sys.stderr)

    cdbg.write(out_file_name,
                cDBGSerialization.enum_from_str(out_format))
    
    if verbose:
        print(f'Finished writing cDBG.', file=sys.stderr)
