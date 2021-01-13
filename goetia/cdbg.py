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

from goetia.serialization import cDBGSerialization


async def compute_connected_component_callback(msg,
                                               cdbg_type,
                                               cdbg,
                                               out_file,
                                               sample_size):
        
    result = cdbg_type.compute_connected_component_metrics(cdbg,
                                                            sample_size)
    n_comps, min_comp, max_comp, size_dist = result

    data = {'t': msg.t,
            'sample_name': msg.sample_name,
            'n_components': n_comps,
            'max': max_comp,
            'min': min_comp,
            'size_dist': size_dist}

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
                                                bins):
    
    counts = cdbg_type.compute_unitig_fragmentation(cdbg,
                                                    bins)
    data = {'t': msg.t,
            'sample_name': msg.sample_name}
    for i in range(len(bins) - 1):
        bin_start, bin_end = bins[i], bins[i+1]
        data[f'[{bin_start},{bin_end})'] = counts[i]
    data[f'[{bin_end},Inf)'] = counts[-1]

    async with curio.aopen(out_file, 'a') as fp:
        if await fp.tell() != 0:
            await fp.write(',\n')
        else:
            await fp.write('[\n')
        await fp.write(json.dumps(data))


async def write_cdbg_metrics_callback(msg,
                                      compactor,
                                      out_file):

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

    async with curio.aopen(out_file, 'a') as fp:
        if await fp.tell() != 0:
            await fp.write(',\n')
        else:
            await fp.write('[\n')
        await fp.write(json.dumps(data))


async def write_cdbg_callback(msg,
                              cdbg,
                              out_file_prefix,
                              out_format):

    out_file_name = f'{out_file_prefix}.{msg.t}.{out_format}'
    cdbg.write(out_file_name,
                cDBGSerialization.enum_from_str(out_format))
