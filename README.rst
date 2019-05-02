.. image:: https://travis-ci.org/camillescott/boink.svg?branch=master
    :target: https://travis-ci.org/camillescott/boink
    
boink
-----

Installation
============

For this repository with anaconda::

    git clone https://github.com/camillescott/boink && cd boink
    git submodule update --init --recursive

    conda create -y -n libboink -c conda-forge python=3 cmake cxx-compiler c-compiler clangdev libcxx libstdcxx-ng libgcc-ng pytest numpy scipy openmp
    conda activate libboink
    pip install -r requirements.txt

    mkdir build; cd build
    cmake ..
    make
    make install

