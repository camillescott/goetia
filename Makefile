LIB_BUILD_DIR = cmake-build
CONDA_FRONTEND = mamba
BUILD_TYPE = Release

# $(LIB_BUILD_DIR)/libgoetiaCppyy.so $(LIB_BUILD_DIR)/libgoetia.so

all: version install-lib 

clean: FORCE
	rm -rf build dist
	rm -rf goetia.egg-info
	rm -rf $(LIB_BUILD_DIR)
	rm -f goetia/libgoetiaCppyy.rootmap
	rm -f goetia/libgoetiaCppyy.so
	rm -f goetia/libgoetiaCppyy_rdict.pcm
	rm -f goetia/goetia.map
	@find ./ -type d -name __pycache__ -exec rm -rf {} +

version: CMAKE_VERSION:= $(shell python version.py --cmake)
version: VERSION:= $(shell python version.py)
version: FORCE
	echo ${CMAKE_VERSION} > include/goetia/VERSION
	echo ${VERSION} > goetia/VERSION 

create-dev-env: environment_dev.yml
	${CONDA_FRONTEND} create -n goetia-dev python=3.8
	${CONDA_FRONTEND} env update -f environment_dev.yml

configure: FORCE
	cmake -H. -B$(LIB_BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -G Ninja

build-lib: configure
	cd $(LIB_BUILD_DIR) && ninja -v

install-lib: build-lib
	cmake --install $(LIB_BUILD_DIR)

bdist_wheel: install-lib
	python setup.py bdist_wheel

install: bdist_wheel
	python -m pip install --no-deps --force-reinstall dist/goetia-`python -c "import sys; sys.stdout.write(open('goetia/VERSION').read().strip()[1:])"`-py3-none-any.whl

dev-install: install-lib
	python -m pip install --no-deps -e .

compile-commands: FORCE
	cmake -H. -Bcmake-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=YES
	ln -fs cmake-debug/compile_commands.json .

test: FORCE
	pytest -v --benchmark-disable tests/

FORCE:
