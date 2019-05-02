.. image:: https://travis-ci.org/camillescott/boink.svg?branch=master
    :target: https://travis-ci.org/camillescott/boink
    
boink
-----

Installation
============

For this repository with anaconda::

    conda create -n libboink python=3 cmake cxx-compiler c-compiler clangdev libcxx libstdcxx-ng libgcc-ng pytest numpy scipy openmp
    conda activate libboink
    pip install -r requirements.txt

    git clone https://github.com/camillescott/boink
    cd boink
    git submodule update --init --recursive
    git checkout cppyy

    mkdir build; cd build
    cmake ..
    make install

