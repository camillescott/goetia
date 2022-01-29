#!/usr/bin/env python
from __future__ import print_function
import argparse
import os, sys, subprocess

import cppyy_backend
from cppyy_backend._get_cppflags import get_cppflags

MYHOME = os.path.dirname(cppyy_backend.__file__)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('header_file', help='Header file to generate reflection info from.')
    parser.add_argument('--extra-cxxflags', default='',
                        help="String of extra compiler flags to pass along "\
                             "to genreflex. MUST come before header_file!")
    parser.add_argument('genreflex_args', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    os.environ['LD_LIBRARY_PATH'] = os.path.join(MYHOME, 'lib')
    genreflex = os.path.join(MYHOME, 'bin', 'genreflex')
    if 'win32' in sys.platform:
        genreflex += '.exe'
    if not os.path.exists(genreflex):
        raise RuntimeError("genreflex not installed in standard location")

    cppyy_extra_flags = get_cppflags() or ''
    extra_flags = ' '.join([cppyy_extra_flags, args.extra_cxxflags])
    print('EXTRA:', extra_flags)
    print('GENREFLEX BIN:', genreflex)
    print('GENREFLEX ARGS:', args.genreflex_args)

    args = [args.header_file, '--cxxflags', extra_flags] + args.genreflex_args
    
    return subprocess.call([genreflex] + args)


if __name__ == "__main__":
    sys.exit(main())
