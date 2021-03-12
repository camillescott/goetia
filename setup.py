#!/usr/bin/env python
import os
import sys

from version import get_version

with open(os.path.join('goetia', 'VERSION')) as fp:
    version = fp.read().strip()

setup(
    version = version
)

