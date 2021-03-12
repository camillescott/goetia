#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : merge-cppyy-maps.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 11.03.2021

import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('map_files', nargs='+')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args()

    # collect maps
    maps = []
    for map_file in args.map_files:
        with open(map_file) as fp:
            payload = json.load(fp)
            maps.append(payload[0])

    # dump map
    json.dump(maps, args.output, indent=1, sort_keys=True)


if __name__ == '__main__':
    main()
