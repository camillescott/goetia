#!/usr/bin/env python

import inspect
import os
import re
import setuptools

from distutils.command.clean import clean
from distutils.util import get_platform
from setuptools.command.build_py import build_py
from setuptools import setup, find_packages
from wheel.bdist_wheel import bdist_wheel


PKG              = 'goetia'
(PKG_NAMESPACE,
 PKG_SIMPLENAME) = PKG.rsplit(".", 1) if '.' in PKG else "", PKG
PKG_DIR          = os.path.dirname(__file__)

class cppyy_build_py(build_py):

    def run(self):
        #
        # Base build.
        #
        build_py.run(self)
        #
        # Custom build.
        #
        #
        # Move CMake output to self.build_lib.
        #
        pkg_subdir = PKG.replace(".", os.path.sep)
        if PKG_NAMESPACE:
            #
            # Implement a pkgutil-style namespace package as per the guidance on
            # https://packaging.python.org/guides/packaging-namespace-packages.
            #
            namespace_init = os.path.join(PKG_NAMESPACE, "__init__.py")
            with open(namespace_init, "w") as f:
                f.write("__path__ = __import__('pkgutil').extend_path(__path__, __name__)\n")
            self.copy_file(namespace_init, os.path.join(self.build_lib, namespace_init))


class cppyy_clean(clean):

    def run(self):
        #
        # Custom clean.
        # TODO: There is no way to reliably clean the "dist" directory.
        #
        #
        #  Base clean.
        #
        clean.run(self)


class cppyy_bdist_wheel(bdist_wheel):

    def finalize_options(self):
        #
        # This is a universal (Python2/Python3), but platform-specific (has
        # compiled parts) package; a combination that wheel does not recognize,
        # thus simply fool it.
        #
        self.plat_name = get_platform()
        bdist_wheel.finalize_options(self)
        self.root_is_pure = True

setup(
    version = '0.15.4.1',
    cmdclass = {
        'build_py': cppyy_build_py,
        'clean':    cppyy_clean
    },

)

