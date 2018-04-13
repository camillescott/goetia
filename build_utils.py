import functools
import os
import shutil
import subprocess
import sys
import sysconfig
import tempfile

from doit.reporter import ConsoleReporter


def get_compiler(CXX):
    output = subprocess.check_output([CXX, '--version'], encoding='UTF-8')
    print(output)
    if 'clang' in output:
        return 'clang'
    if 'GCC' in output:
        return 'gcc'


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


def replace_ext(filename, ext):
    return os.path.splitext(filename)[0] + ext


def distutils_dir_name(dname):
    """Returns the name of a distutils build directory"""
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,
                    platform=sysconfig.get_platform(),
                    version=sys.version_info)


def build_dir():
    return os.path.join("build", distutils_dir_name("temp"))


def lib_dir():
    return os.path.join("build", distutils_dir_name("lib"))


def title_with_actions(task):
    """return task name task actions"""
    if task.actions:
        title = "\n\t".join([str(action) for action in task.actions])
    # A task that contains no actions at all
    # is used as group task
    else:
        title = "Group: %s" % ", ".join(task.task_dep)
    return "%s\n%s"% (task.name, title)


class BoinkReporter(ConsoleReporter):

    def divider(self):
        self.write('\n')
        self.write('=' * 80)
        self.write('\n')

    def execute_task(self, task):
        if task.actions and (task.name[0] != '_'):
            self.divider()
            super().execute_task(task)

    def skip_uptodate(self, task):
        if task.name[0] != '_':
            self.divider()
            super().skip_uptodate(task)
    
    def skip_ignore(self, task):
        self.divider()
        super().skip_ignore(task)


def get_cpp_includes(filename):
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('#include'):
                line = line.split()
                name_token = line[1]
                include_filename = name_token.strip('"<>"')
                if include_filename.endswith('.hh'):
                    yield include_filename


def get_pxd_includes(filename):
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('cdef extern from'):
                tokens = line.split()
                include_filename = tokens[3].strip('"')
                include_filename = os.path.basename(include_filename)
                if include_filename.endswith('.hh'):
                    yield include_filename


def get_object_dependencies(source_filename):
    header_filename = os.path.splitext(source_filename)[0] + '.hh'
    for fn in get_cpp_includes([source_filename, header_filename]):
        yield fn


def get_pyx_imports(filename, mods, MOD_EXT):
    mods = [os.path.basename(fn) for fn in mods]
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if 'cimport'  in line:
                tokens = line.split()
                mod = tokens[1].strip('.')
                mod = mod.split('.')[-1]
                mod = mod + MOD_EXT
                if mod in mods:
                    yield mod


def resolve_dependencies(PYX_FILES, PYX_NAMES, PXD_FILES, PXD_NAMES,
                         SOURCE_FILES, SOURCE_NAMES, HEADER_FILES, HEADER_NAMES,
                         EXT_SOS, MOD_EXT, INCLUDE_DIR, PKG):
    HEADERS = dict(zip(HEADER_NAMES, HEADER_FILES))
    PXD     = dict(zip(PXD_NAMES, PXD_FILES))
    PYX     = dict(zip(PYX_NAMES, PYX_FILES))
    SOURCES = dict(zip(SOURCE_NAMES, SOURCE_FILES))

    def _resolve_dependencies(filename, includes):
        if not os.path.isfile(filename):
            print('not found:', filename)
            return includes

        if filename.endswith('.hh'):
            _includes = set()
            _includes.update(get_cpp_includes(filename))
            _includes = {HEADERS[inc] for inc in _includes}
            return includes | _includes
        elif filename.endswith('.cc'):
            _includes = set()
            _header_path = HEADERS[replace_ext(os.path.basename(filename), '.hh')]
            _includes.add(_header_path)
            _includes.update(_resolve_dependencies(_header_path, _includes))
            _includes.update(get_cpp_includes(filename))
            return includes | _includes
        elif filename.endswith('.pyx'):
            _includes = set()
            _mod_imports = get_pyx_imports(filename, EXT_SOS, MOD_EXT)
            _mod_imports = [os.path.join(PKG, fn) for fn in _mod_imports]
            _includes.update(_mod_imports)
            _pxd_path = replace_ext(filename, '.pxd')
            if os.path.isfile(_pxd_path):
                _includes.add(_pxd_path)
            _includes.update(_resolve_dependencies(_pxd_path, _includes))
            return includes | _includes
        elif filename.endswith('.pxd'):
            _includes = set()
            _mod_imports = get_pyx_imports(filename, EXT_SOS, MOD_EXT)
            _mod_imports = [os.path.join(PKG, fn) for fn in _mod_imports]
            _includes.update(_mod_imports)
            _headers = [HEADERS[h] for h in get_pxd_includes(filename)]
            _includes.update(_headers)
            return includes | _includes
        else:
            return includes

    DEP_MAP = {}
    for mod in EXT_SOS:
        pyx = os.path.basename(mod).split('.')[0] + '.pyx'
        cpp = replace_ext(pyx, '.cpp')
        os.path.join(PKG, cpp)
        pyx_file = PYX[pyx]
        includes = {pyx_file}
        includes.update(_resolve_dependencies(pyx_file, includes))
        DEP_MAP[mod] = includes
        for fn, deps in DEP_MAP.items():
            if fn in deps:
                deps.remove(fn)

    return DEP_MAP
