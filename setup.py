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
from setuptools import Extension
from distutils.spawn import spawn
from distutils.sysconfig import get_config_vars
from distutils.dist import Distribution
from distutils.errors import DistutilsPlatformError

import numpy

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext as _build_ext
except ImportError:
    print('Cython not found.', file=sys.stderr)
    from setuptools.command.build_ext import build_ext as _build_ext

# strip out -Wstrict-prototypes; a hack suggested by
# http://stackoverflow.com/a/9740721
# proper fix coming in http://bugs.python.org/issue1222585
# numpy has a "nicer" fix:
# https://github.com/numpy/numpy/blob/master/numpy/distutils/ccompiler.py
OPT = get_config_vars('OPT')[0]
os.environ['OPT'] = " ".join(
    flag for flag in OPT.split() if flag != '-Wstrict-prototypes'
)

# Checking for OpenMP support. Currently clang doesn't work with OpenMP,
# so it needs to be disabled for now.
# This function comes from the yt project:
# https://bitbucket.org/yt_analysis/yt/src/f7c75759e0395861b52d16921d8ce3ad6e36f89f/yt/utilities/lib/setup.py?at=yt


def check_for_openmp():
    """Check for OpenMP support."""
    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    exit_code = 1

    if os.name == 'nt':
        return False

    try:
        os.chdir(tmpdir)

        # Get compiler invocation
        compiler = os.getenv('CC', 'cc')

        # Attempt to compile a test script.
        # See http://openmp.org/wp/openmp-compilers/
        filename = r'test.c'
        source = open(filename, 'wt', 1)
        source.write(
            """
            #include <omp.h>
            #include <stdio.h>
            int main() {
            #pragma omp parallel
            printf("Hello from thread %d, nthreads %d",
                    omp_get_thread_num(), omp_get_num_threads());
            }
            """
        )
        with open(os.devnull, 'w') as fnull:
            exit_code = subprocess.call([compiler, '-fopenmp', filename],
                                        stdout=fnull, stderr=fnull)

        # Clean up
        source.close()
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    return exit_code == 0

def distutils_dir_name(dname):
    """Returns the name of a distutils build directory"""
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,
                    platform=sysconfig.get_platform(),
                    version=sys.version_info)

def build_dir():
    return path_join("build", distutils_dir_name("temp"))


# Don't forget to update lib/Makefile with these flags!
EXTRA_COMPILE_ARGS = ['-O3', '-std=c++14', '-pedantic']
EXTRA_LINK_ARGS = ['--verbose']

if sys.platform == 'darwin':
    # force 64bit only builds
    EXTRA_COMPILE_ARGS.extend(['-arch', 'x86_64', '-mmacosx-version-min=10.7',
                               '-stdlib=libc++'])
else:
    EXTRA_COMPILE_ARGS.append('-fdiagnostics-color')

if check_for_openmp():
    EXTRA_COMPILE_ARGS.extend(['-fopenmp'])
    EXTRA_LINK_ARGS.extend(['-fopenmp'])

EXTENSION_MODS = []
for cython_ext in glob.glob(os.path.join("boink", "*.pyx")):

    CY_EXTENSION_MOD_DICT = \
        {
            "sources": [cython_ext],
            "extra_compile_args": EXTRA_COMPILE_ARGS,
            "extra_link_args": EXTRA_LINK_ARGS,
            "depends": [],
            "libraries": ['oxli'],
            "include_dirs": [numpy.get_include()],
            "language": "c++"
        }
    
    ext_name = "boink.{0}".format(splitext(os.path.basename(cython_ext))[0])
    CY_EXTENSION_MOD = Extension(ext_name, ** CY_EXTENSION_MOD_DICT)
    EXTENSION_MODS.append(CY_EXTENSION_MOD)

SCRIPTS = []
SCRIPTS.extend([path_join("scripts", script)
                for script in os_listdir("scripts")
                if script.endswith(".py")])

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.3",
    "Programming Language :: Python :: 3.4",
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
        "packages": ['boink'],
        "package_data": {'boink/': ['*.pxd']},
        "install_requires": ['screed >= 0.9', 'bz2file'],
        "setup_requires": ["pytest-runner>=2.0,<3dev",
                           'Cython>=0.25.2', "setuptools>=18.0"],
        "scripts": SCRIPTS,
        # "entry_points": { # Not ready for distribution yet.
        #    'console_scripts': [
        #        "oxli = oxli:main"
        #    ]
        # },
        "ext_modules": cythonize(EXTENSION_MODS, gdb_debug = True),
        #                        compiler_directives = {"language_level": 3}),
        # "platforms": '', # empty as is conveyed by the classifiers below
        # "license": '', # empty as is conveyed by the classifier below
        "include_package_data": True,
        "zip_safe": False,
        "classifiers": CLASSIFIERS
    }



setup(**SETUP_METADATA)
