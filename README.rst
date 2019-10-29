.. image:: https://travis-ci.org/camillescott/boink.svg?branch=master
    :target: https://travis-ci.org/camillescott/boink

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/camillescott/boink/master?filepath=examples%2FStreaming%20Sourmash%20Demo.ipynb
    
boink
-----

Installation
============

Conda
~~~~~

We recommend using `conda <https://docs.conda.io/en/latest/miniconda.html>`_. Within a conda
environment, install with::

    conda install graph-boink

This will install the boink python package, and install the `libboink` shared library
and its headers into your conda prefix.

Building from Source
~~~~~~~~~~~~~~~~~~~~

To build and install from source, first clone the repo::

    git clone https://github.com/camillescott/boink && cd boink
    git submodule update --init --recursive

Then create a conda environment::

    conda create -y -n libboink -c conda-forge python=3 cppyy cmake cxx-compiler c-compiler clangdev libcxx libstdcxx-ng libgcc-ng pytest numpy scipy openmp python-clang screed blessings pytest-benchmark pyfiglet py-cpuinfo
    conda activate libboink

Then build and install with cmake::

    mkdir build; cd build
    cmake ..
    make
    make install
