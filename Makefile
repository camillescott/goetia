LIB_BUILD_DIR = cmake_build

# $(LIB_BUILD_DIR)/libgoetiaCppyy.so $(LIB_BUILD_DIR)/libgoetia.so

all: version install-lib 

clean: FORCE
	rm -rf build dist
	rm -rf goetia.egg-info
	rm -rf $(LIB_BUILD_DIR)
	@find ./ -type d -name __pycache__ -exec rm -rf {} +

version: FORCE
	VERSION=`python version.py`
	CMAKE_VERSION=`python version.py --cmake`
	echo $(VERSION) > goeta/VERSION
	echo $(CMAKE_VERSION) > include/goetia/VERSION

create-dev-env: environment_dev.yml
	mamba create -n goetia-dev python=3.8
	mamba env update -f environment_dev.yml

build-lib: FORCE
	cmake -H. -B$(LIB_BUILD_DIR) -G Ninja
	cmake --build $(LIB_BUILD_DIR) -- -v

install-lib: build-lib
	cmake --install $(LIB_BUILD_DIR)

bdist_wheel: install-lib
	python setup.py bdist_wheel

install: bdist_wheel
	python -m pip install --no-deps --force-reinstall dist/goetia-`python -c "import sys; sys.stdout.write(open('goetia/VERSION').read().strip()[1:])"`-py3-none-any.whl

dev-install: install-lib
	python -m pip install --no-deps -e .

compile-commands: FORCE
	cmake -H. -Bcmake_debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=YES
	ln -fs cmake_debug/compile_commands.json .

FORCE:
