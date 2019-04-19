For this repository with anaconda::

    conda create -n libboink python=3 cmake cxx-compiler c-compiler clangdev libcxx libstdcxx-ng libgcc-ng pytest numpy scipy
    conda activate libboink
    pip install -r requirements.txt

    git clone https://github.com/camillescott/boink
    cd boink
    git submodule update --init --recursive
    git checkout cppyy

    mkdir build; cd build
    cmake ..
    make install

