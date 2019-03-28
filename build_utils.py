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

def walk_ext(root, ext='.cc'):
    result = []
    for root, dirs, files in os.walk(root):
        for fn in files:
            if fn.endswith(ext):
                result.append((root, fn))
    return sorted(result, key=lambda p: (os.path.basename(p[0]), replace_ext(p[1], '')))

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
        title += "\nDeps:    " + ", ".join([str(dep) for dep in task.file_dep])
        title += "\nTargets: " + ", ".join([str(tgt) for tgt in task.targets])
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
        self.write(result['task'].title())
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


def get_gcc_includes(source, include_directives, pkg, ignore=[]):
    from sh import gcc
    includes =  flatten([ln.strip('\ ').split() for ln in \
                        gcc(*include_directives, '-std=c++14', '-MM', '-MG', source).split('\n')])[1:]
    if not ignore:
        return includes

    result = []
    for include in includes:
        skip = False
        for path in ignore:
            if include.startswith(path):
                skip = True
        if skip is False:
            result.append(include)
    return result


