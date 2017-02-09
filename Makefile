#CMD=g++ -shared -L/home/camille/miniconda/envs/dev.py3.khmer.assembly/lib 
#-Wl,-rpath=/home/camille/miniconda/envs/dev.py3.khmer.assembly/lib 
#-lpython3.5m -L/home/camille/miniconda/envs/dev.py3.khmer.assembly/lib 
#-I/home/camille/miniconda/envs/dev.py3.khmer.assembly/include/python3.5m 
#-I/work/khmer/khmer -I/work/khmer/lib -std=c++11 -fPIC test.c

PKG=boink
MODEXT=$(shell python -c \
       "import sysconfig;print(sysconfig.get_config_var('SO'))")

MODS=$(PKG)/*$(MODEXT)

all:
	python setup.py build_ext --inplace



clean: FORCE
	rm -f $(PKG)/*.cpp
	@find ./ -type d -name __pycache__ -exec rm -rf {} +
	@find ./$(PKG)/ -type f -name *$(MODEXT) -exec rm -f {} +
	@find ./$(PKG)/ -type f -name *.pyc -exec rm -f {} +
	@find ./$(PKG)/ -type f -name *.cpp -exec rm -f {} +
	rm -rf build dist $(PKG).egg-info

FORCE:


# Use this to print the value of a Makefile variable
# Example `make print-VERSION`
# From https://www.cmcrossroads.com/article/printing-value-makefile-variable
print-%  : ; @echo $* = $($*)
