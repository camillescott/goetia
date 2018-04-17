#! /usr/bin/env python
# vim: set fileencoding=utf-8

import glob
import os
import sys
from os import listdir as os_listdir
from os.path import join as path_join
from os.path import splitext
import shutil
import subprocess
import sys
import sysconfig
import tempfile

from setuptools import setup
from distutils.spawn import spawn
from distutils.sysconfig import get_config_vars
from distutils.dist import Distribution
from distutils.errors import DistutilsPlatformError
from setuptools.command.build_ext import build_ext as _build_ext

import numpy
from build_utils import (check_for_openmp,
                         distutils_dir_name,
                         build_dir,
                         get_compiler)


SCRIPTS = glob.glob(os.path.join('scripts', '*'))
MOD_EXT = sysconfig.get_config_var('SO')
MODS    = glob.glob(os.path.join('boink', '*'+MOD_EXT))

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]


SETUP_METADATA = \
    {
        "name": "boink",
        "version": "0.1",
        "description": '',
        "long_description": open("README.rst").read(),
        "author": "Camille Scott",
        "author_email": 'camille.scott.w@gmail.com',
        "packages": ['boink', 'boink.tests'],
        "package_data": {'boink/': ['*.pxd', '*.pxi'] + MODS},
        "install_requires": ['screed >= 0.9', 'bz2file'],
        "setup_requires": ["pytest-runner>=2.0,<3dev",
                           'Cython>=0.25.2', "setuptools>=18.0"],
        "scripts": SCRIPTS,
        "include_package_data": True,
        "zip_safe": False,
        "classifiers": CLASSIFIERS
    }

class BoinkBuildExt(_build_ext):

    def run(self):
        '''
        if HAS_CYTHON:
            print('*** NOTE: Found Cython, extension files will be '
                  'transpiled if this is an install invocation.')
        else:
            print('*** WARNING: Cython not found, assuming cythonized '
                  'files available for compilation.')

        extensions = ('{0}:{1}'.format(x, y) for x, y in EXTENSION_NAMES)
        print('*** EXTENSIONS:', ', '.join(extensions))
        print('*** INCLUDES:', ', '.join(DEPENDS))
        print('*** SOURCES:', ', '.join(SOURCES))

        _build_ext.run(self)
        '''
        _build_ext.run(self)

setup(**SETUP_METADATA,
      cmdclass = {
         'build_ext': BoinkBuildExt
      })
