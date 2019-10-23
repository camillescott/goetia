# This program is placed into the public domain.

"""
Gets the current version number.
If in a git repository, it is the current git tag.
Otherwise it is the one contained in the PKG-INFO file.
To use this script, simply import it in your setup.py file
and use the results of get_version() as your package version:
    from version import *
    setup(
        ...
        version=get_version(),
        ...
    )
"""

__all__ = ('get_version')

import argparse
from os.path import dirname, isdir, join
import os
import re
import subprocess

version_re = re.compile('^Version: (.+)$', re.M)


def get_version(cmake=False):
    d = dirname(__file__)

    version = []

    if isdir(join(d, '.git')):
        # Get the version using "git describe".
        cmd = 'git describe --tags'.split()
        try:
            tag = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get version number from git tags')
            exit(1)
        tag = tag.strip('v')

        vlevels = tag.split('.')

        # PEP 386 compatibility
        if '-' in vlevels[-1]:
            postname = 'post'.join(tag[-1].split('-')[:2])
            version.extend(vlevels[:-1])
            versiona.append(postname)
        else:
            version = vlevels

        # Don't declare a version "dirty" merely because a time stamp has
        # changed. If it is dirty, append a ".dev1" suffix to indicate a
        # development revision after the release.
        with open(os.devnull, 'w') as fd_devnull:
            subprocess.call(['git', 'status'],
                            stdout=fd_devnull, stderr=fd_devnull)

        cmd = 'git diff-index --name-only HEAD'.split()
        try:
            dirty = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get git index status')
            exit(1)

        if dirty != '':
            version.append('dev1')

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)


    print(version)
    return '.'.join(version)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmake', action='store_true', default=False)
    args = parser.parse_args()

    print(get_version(cmake=args.cmake))
