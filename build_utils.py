import functools
import itertools
import os
import shutil
from shutil import rmtree
import subprocess
import sys
import sysconfig
import tempfile

from doit.reporter import ConsoleReporter
from jinja2 import Environment, PackageLoader, select_autoescape

'''
*
* Compiler stuff
*
'''
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

'''
*
* File stuff
*
'''

def flatten(L):
    return [item for sublist in L for item in sublist]


def replace_ext(filename, ext):
    return os.path.splitext(filename)[0] + ext


def replace_exts(files, ext):
    if type(files) is list:
        return [replace_ext(fn, ext) for fn in files]
    else:
        return {k:replace_ext(v, ext) for k,v in files.items()}


def isfile_filter(files):
    if type(files) is list:
        return [fn for fn in files if os.path.isfile(fn)]
    else:
        return {k:v for k,v in files.items() if os.path.isfile(v)}



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

def clean_folder(target):
    '''Function for doit task's `clean` parameter to remove a folder.

    Args:
        target (str): The folder to remove.
    '''
    print('Remove folder', target)
    try:
        rmtree(target)
    except OSError:
        pass

'''
*
* Printing stuff
*
'''

class TermCodes:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def text(txt, code):
    return code + txt + TermCodes.ENDC

def get_console_size():
    try:
        rows, columns = os.popen('stty size', 'r').read().split()
    except:
        rows, columns = 80, 80
    return int(rows), int(columns)


def major_divider(text=''):
    _, columns = get_console_size()
    if text:
        text = '== {0} '.format(text)
    return '\n{0}{1}'.format(text, '=' * (columns-len(text)))


def minor_divider(text=''):
    _, columns = get_console_size()
    if text:
        text = '-- {0} '.format(text)
    return '{0}{1}'.format(text, '-' * (columns-len(text)))


def title_with_actions(task):
    """return task name task actions"""
    if task.actions:
        title = "\n".join([str(action) for action in task.actions])
        title += "\nDeps: " + ", ".join([str(dep) for dep in task.file_dep])
    # A task that contains no actions at all
    # is used as group task
    else:
        title = "Group: %s" % ", ".join(task.task_dep)
    return "Name: {name}\n{title}\n".format(
            name=text(task.name, TermCodes.HEADER),
            title=title)


class BoinkReporter(ConsoleReporter):

    def execute_task(self, task):
        if task.actions and (task.name[0] != '_'):
            self.write(text(major_divider('RUN TASK'), TermCodes.WARNING))
            self.write('\n')
            self.write(task.title())
            self.write(text(minor_divider('task output'),
                       TermCodes.OKGREEN))
            self.write('\n')

    def skip_uptodate(self, task):
        if task.name[0] != '_':
            self.write(text(major_divider('TASK UP-TO-DATE'),
                            TermCodes.OKBLUE))
            self.write('\n')
            self.write(task.title())
    
    def skip_ignore(self, task):
        self.write(major_divider())
        super().skip_ignore(task)

    def _write_failure(self, result, write_exception=True):
        self.write(text(major_divider('TASK FAILURE'), TermCodes.FAIL))
        err = text(result['exception'].get_name(), TermCodes.FAIL)
        tsk = text(result['task'].name, TermCodes.HEADER)
        msg = '%s - taskid:%s\n' % (err, tsk)
        self.write(msg)
        if write_exception:
            self.write(result['exception'].get_msg())
            self.write("\n")

'''
*
* Cython stuff
*
'''


def get_gcc_includes(source, includes, pkg):
    from sh import gcc
    res = gcc(*includes, '-std=c++14', '-M', source)
    return [os.path.abspath(token) for token in res.split() if pkg in token]


def get_cpp_includes(filename, include_dir):
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('#include'):
                line = line.split()
                name_token = line[1]
                include_filename = name_token.strip('"<>"')
                if include_filename.endswith('.hh'):
                    yield os.path.join(include_dir, include_filename)


def get_pxd_includes(filename, include_dir):
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('cdef extern from'):
                tokens = line.split()
                include_filename = tokens[3].strip('"')
                include_filename = os.path.basename(include_filename)
                include_filename = os.path.join(include_dir, include_filename)
                if include_filename.endswith('.hh') and \
                   os.path.isfile(include_filename):
                    yield include_filename


def get_object_dependencies(source_filename):
    header_filename = os.path.splitext(source_filename)[0] + '.hh'
    for fn in get_cpp_includes([source_filename, header_filename]):
        yield fn


def get_pyx_imports(filename, MOD_FILES, PKG):
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if 'cimport'  in line:
                tokens = line.split()
                mod = tokens[1].strip('.')
                mod = mod.split('.')[-1]
                if mod in MOD_FILES:
                    pyx = os.path.join(PKG, mod+'.pyx')
                    pxd = replace_ext(pyx, '.pxd')
                    yield from  isfile_filter([pyx, pxd])


def resolve_dependencies(PYX_FILES, PXD_FILES,
                         MOD_FILES, INCLUDE_DIR, PKG):

    def _resolve_dependencies(filename, includes):
        if not os.path.isfile(filename):
            return includes

        if filename.endswith('.hh'):
            _includes = set()
            _includes.update(get_cpp_includes(filename, INCLUDE_DIR))
            return includes | _includes
        elif filename.endswith('.cc'):
            _includes = set()
            _header_path = replace_ext(os.path.basename(filename), '.hh')
            _includes.add(_header_path)
            _includes.update(_resolve_dependencies(_header_path, _includes))
            _includes.update(get_cpp_includes(filename, INCLUDE_DIR))
            return includes | _includes
        elif filename.endswith('.pyx'):
            _includes = set()
            _mod_imports = get_pyx_imports(filename, MOD_FILES, PKG)
            _includes.update(_mod_imports)
            _pxd_path = replace_ext(filename, '.pxd')
            if os.path.isfile(_pxd_path):
                _includes.add(_pxd_path)
            _includes.update(_resolve_dependencies(_pxd_path, _includes))
            return includes | _includes
        elif filename.endswith('.pxd'):
            _includes = set()
            _mod_imports = get_pyx_imports(filename, MOD_FILES, PKG)
            _includes.update(_mod_imports)
            _includes.update(get_pxd_includes(filename, INCLUDE_DIR))
            return includes | _includes
        else:
            return includes

    DEP_MAP = {}
    for mod, soname in MOD_FILES.items():
        pyx_file = PYX_FILES[mod]
        cpp = replace_ext(pyx_file, '.cpp')
        includes = {pyx_file}
        includes.update(_resolve_dependencies(pyx_file, includes))
        DEP_MAP[mod] = includes
        #for fn, deps in DEP_MAP.items():
        #    if fn in deps:
        #        deps.remove(fn)
    for mod, mod_deps in DEP_MAP.items():
        mod_path = os.path.join(PKG, MOD_FILES[mod])
        if mod_path in mod_deps:
            mod_deps.remove(mod_path)
    return DEP_MAP

'''
*
* Jinja stuff
*
'''

def get_template_env():
    return Environment(loader=PackageLoader('boink', 'templates'),
                       trim_blocks=True,
                       lstrip_blocks=True)


def generate_cpp_params(type_dict, cppclass='class'):
    results = []
    cppclasses = []
    for param_bundle in itertools.product(*(v for _, v in type_dict.items())):
        result = {}
        for name, param in zip(type_dict.keys(), param_bundle):
            result[name] = param
        result['suffix'] = '_'.join(param_bundle)
        result['params'] = ','.join(param_bundle)
        results.append(result)
        cppclasses.append('{0}[{1}]'.format(cppclass, result['params']))
    return results, cppclasses


def get_templates(mod_prefix):
    env = get_template_env()

    try:
        pxd_template = env.get_template('{0}.tpl.pxd'.format(mod_prefix))
    except Exception as e:
        raise RuntimeError("Error parsing {0}: {1}".format('{0}.tpl.pxd'.format(mod_prefix),
                                                           e))
    pxd_dst_path = os.path.join('boink', pxd_template.name + '.pxi')

    try:
        pyx_template = env.get_template('{0}.tpl.pyx'.format(mod_prefix))
    except Exception as e:
        raise RuntimeError("Error parsing {0}: {1}".format('{0}.tpl.pyx'.format(mod_prefix),
                                                           e))
    pyx_dst_path = os.path.join('boink', pyx_template.name + '.pxi')

    return (pxd_template, pxd_dst_path,
            pyx_template, pyx_dst_path)


def render(pxd_template, pxd_dst_path,
           pyx_template, pyx_dst_path, type_bundles):
    
    with open(pxd_dst_path, 'w') as fp:
        try:
            rendered = pxd_template.render(dst_filename=pxd_dst_path,
                                           tpl_filename=pxd_template.name,
                                           type_bundles=type_bundles)
            fp.write(rendered)
        except Exception as e:
            raise RuntimeError("Error rendering {0}: {1}".format(pxd_template.name,
                                                                 e))
    with open(pyx_dst_path, 'w') as fp:
        try:
            rendered = pyx_template.render(dst_filename=pyx_dst_path,
                                           tpl_filename=pyx_template.name,
                                           type_bundles=type_bundles)
            fp.write(rendered)
        except Exception as e:
            raise RuntimeError("Error rendering {0}: {1}".format(pyx_template.name,
                                                                 e))
 
