#!/usr/bin/env doit

from collections import OrderedDict
import glob
import itertools
import os
from pprint import pprint
import sysconfig
import sys

from doit import create_after, get_var
from doit.task import clean_targets
from doit.tools import run_once
import yaml

from build_utils import (check_for_openmp,
                         clean_folder,
                         distutils_dir_name,
                         build_dir,
                         lib_dir,
                         BoinkReporter,
			 title_with_actions,
                         replace_ext,
                         replace_exts,
                         isfile_filter,
                         flatten,
                         get_gcc_includes,
                         walk_ext)

DOIT_CONFIG  = {'verbosity': 2,
                'reporter': BoinkReporter,
                'default_tasks': ['libboink',
                                  ]}

#
# General package information
#
PKG          = 'boink'
VERSION      = open(os.path.join(PKG, 'VERSION')).read().strip()

MOD_EXT      = sysconfig.get_config_var('SO')
PY_FILES     = list(glob.glob('boink/**/*.py', recursive=True))

#
# libboink vars
#
INCLUDE_DIR  = os.path.join('include', PKG)
SOURCE_DIR   = os.path.join('src', PKG)

HEADER_FILES = walk_ext(INCLUDE_DIR, '.hh')
HEADER_NAMES = [header for d, header in HEADER_FILES]

SOURCE_FILES = [(d, f) for d, f in walk_ext(SOURCE_DIR, '.cc') \
                if not d.endswith('benchmarks')]
SOURCE_NAMES = [source for d, source in SOURCE_FILES]

LIB_VERSION = '1'
LIBBUILDDIR = '_libbuild'
#LIBDIR = os.path.join('src', PKG)

CXX         = get_var('CXX', os.environ.get('CXX', 'c++'))
INCLUDES    = ['-I' + os.path.abspath('include/'), '-I.', '-I' + os.path.abspath('include/sourmash/')]

# for anaconda
PREFIX      = get_var('PREFIX', sysconfig.get_config_var('prefix'))
if PREFIX is not None:
    INCLUDES += ['-I'+os.path.join(PREFIX, 'include')]

if sys.platform == 'darwin':
    SHARED_EXT   = 'dylib'
    SONAME       = 'lib{0}.{1}.{2}'.format(PKG, SHARED_EXT, LIB_VERSION)
    SONAME_FLAGS = ['-install_name', os.path.join(PREFIX, 'lib', SONAME),
                    '-compatibility_version', LIB_VERSION,
                    '-current_version', LIB_VERSION]
else:
    SHARED_EXT   = 'so'
    SONAME       = 'lib{0}.{1}.{2}'.format(PKG, SHARED_EXT, LIB_VERSION)
    SONAME_FLAGS = ['-Wl,-soname={0}'.format(SONAME)]

LIBBOINKSO       = os.path.join(LIBBUILDDIR, SONAME)
LIBBOINKSHARED   = os.path.join('lib', 'lib{0}.{1}'.format(PKG, SHARED_EXT))
LIBBOINKBUILDDIR = os.path.join(LIBBUILDDIR, PKG)

OBJECT_FILES = [(os.path.join(LIBBOINKBUILDDIR, d.replace(SOURCE_DIR, '').strip('/')),
                 replace_ext(source, '.o')) \
                for d, source in SOURCE_FILES]

# Start libboink compile vars

PROFILING    = get_var('PROFILING', False)
PROFILER     = 'gprof'
COLOR        = True

DEBUG_CDBG   = get_var('DEBUG_CDBG', False)
DEBUG_CPTR   = get_var('DEBUG_CPTR', False)
DEBUG_ASMLY  = get_var('DEBUG_ASMY', False)
DEBUG_EVENTS = get_var('DEBUG_EVENTS', False)
DEBUG_ALL    = get_var('DEBUG_ALL', False)
DEBUG        = get_var('DEBUG', DEBUG_CDBG or \
                                DEBUG_CPTR or \
                                DEBUG_ALL or \
                                DEBUG_ASMLY or \
                                DEBUG_EVENTS)

HEAP_PROFILE = get_var('HEAP_PROFILE', False)

if DEBUG_ALL:
    DEBUG        = True
    DEBUG_CDBG   = True
    DEBUG_CPTR   = True
    DEBUG_ASMLY  = True
    DEBUG_EVENTS = True

DEBUG_FLAGS = ['-gdwarf']

WARNINGS    = ['-Wall']
COMMON      = ['-O3', '-fPIC', '-fno-omit-frame-pointer']

CPPFLAGS    = []

CXXFLAGS    = COMMON + WARNINGS
CXXFLAGS   += ['-Wstrict-null-sentinel',
               '-std=c++14',
               '-frecord-gcc-switches']
CXXFLAGS   += INCLUDES
CXXFLAGS   += CPPFLAGS

CFLAGS      = ['-Wshadow', '-Wcast-align', '-Wstrict-prototypes']
CFLAGS     += INCLUDES
CFLAGS     += CPPFLAGS

LDFLAGS     = []
if PREFIX is not None:
    LDFLAGS += ['-Wl,-rpath,{0}'.format(os.path.join(PREFIX, 'lib')),
                '-L'+os.path.join(PREFIX, 'lib')]
LDFLAGS     += ['-Wl,-rpath,{0}'.format(os.path.abspath('lib')),
                '-Llib', '-lgfakluge', '-lprometheus-cpp-pull',
                '-lprometheus-cpp-push', '-lprometheus-cpp-core',
                '-lpthread', '-lz']


if DEBUG:
    CXXFLAGS += DEBUG_FLAGS
    CFLAGS   += DEBUG_FLAGS
    try:
        CXXFLAGS.remove('-O3')
        CXXFLAGS.append('-O0')
        CXXFLAGS.append('-fkeep-inline-functions')
    except ValueError:
        pass
    if DEBUG_CPTR:
        CXXFLAGS += ['-DDEBUG_CPTR']
    if DEBUG_CDBG:
        CXXFLAGS += ['-DDEBUG_CDBG']
    if DEBUG_ASMLY:
        CXXFLAGS += ['-DDEBUG_ASMLY']
    if DEBUG_EVENTS:
        CXXFLAGS += ['-DDEBUG_EVENTS']

if PROFILING:
    if PROFILER == 'gprof':
        CXXFLAGS += ['-pg']
        CFLAGS   += ['-pg']
        LDFLAGS  += ['-pg']

if HEAP_PROFILE:
    LDFLAGS  += ['-ltcmalloc']

if sys.platform == 'linux':
    LDFLAGS  += ['-pthread']

if check_for_openmp():
    LDFLAGS  += ['-fopenmp']
    CXXFLAGS += ['-fopenmp']

if COLOR:
    CXXFLAGS += ['-fdiagnostics-color']

if sys.platform == 'darwin':
    CXXFLAGS += [ '-mmacosx-version-min=10.7' ]

CXXFLAGS  += ['-DVERSION={0}'.format(VERSION)]



#
# Compilation commands
#
def cxx_command(source, target):
    cmd = [ CXX ] \
          + LDFLAGS \
          + CXXFLAGS \
          + ['-c -o', target, source]
    return ' '.join(cmd)


def link_command(objects, target):
    cmd = [ CXX ] \
          + LDFLAGS \
          + CXXFLAGS \
          + ['-shared', '-o', target] \
          + objects
    return ' '.join(cmd)



#########
#
# doit tasks.
#
#########

def task_display_libboink_config():

    def print_config():
        print('CXX:', CXX)
        print('SOURCES:')
        pprint(list((os.path.join(*source) for source in SOURCE_FILES)))
        print('HEADERS:')
        pprint(list((os.path.join(*header) for header in HEADER_FILES)))
        print('OBJECTS:')
        pprint(list((os.path.join(*obj)    for obj    in OBJECT_FILES)))
        print('libboink:', LIBBOINKSO, LIBBOINKSHARED)
        print('Compile:', cxx_command('SOURCE.cc', 'TARGET.o'))
        print('Link:', link_command(['OBJECT_1.o', 'OBJECT_2.o',
                                     'OBJECT_N.o'], SONAME))
        print('PY_LDSHARED:', sysconfig.get_config_var('LDSHARED'))
    return {'actions':  [print_config],
            'uptodate': [False]}

#
# Tasks preparing for libboink build
#
def task_create_libboink_build_dirs():
    dirs  = set((directory for directory, _ in OBJECT_FILES))
    files = [os.path.join(d, 'touched') for d in dirs]
    
    return {'title':    title_with_actions,
            'actions':  ['mkdir -p {0}'.format(' '.join(dirs)),
                         'touch {0}'.format(' '.join(files))],
            'task_dep': ['display_libboink_config'],
            'targets':  files,
            'clean':    [clean_targets],
            'uptodate': [run_once]}


def task_deploy_prometheus():
    build_cmd   = 'mkdir -p third-party/prometheus-cpp/_build && cd third-party/prometheus-cpp/_build && '\
                  'cmake -DBUILD_SHARED_LIBS=ON .. && make -j 4 && mkdir -p ../_deploy && make DESTDIR=../_deploy install'
    deploy_root = 'third-party/prometheus-cpp/_deploy/usr/local'
    libs        = ['libprometheus-cpp-core.{}'.format(SHARED_EXT),
                   'libprometheus-cpp-pull.{}'.format(SHARED_EXT),
                   'libprometheus-cpp-push.{}'.format(SHARED_EXT)]
    #libs = [os.path.join('lib', lib) for lib in libs]
    return {'title':    title_with_actions,
            'targets':  [os.path.join('lib', lib) for lib in libs] +
                        [os.path.join(deploy_root, 'lib', lib) for lib in libs],
            'actions':  [build_cmd,
                         'mkdir -p lib/',
                         'cp -r {0} {1}'.format(os.path.join(deploy_root, 'include', 'prometheus'),
                                                'include/'),
                         'cp -r {0}/*.{1} {2}'.format(os.path.join(deploy_root, 'lib*'),
                                                     SHARED_EXT,
                                                     'lib/')],
            'task_dep': ['display_libboink_config'],
            'uptodate': [run_once],
            'clean':    [(clean_folder, ['include/prometheus']),
                         clean_targets]}


def task_deploy_gfakluge():
    build_cmd  = 'cd third-party/gfakluge && make'
    cp_lib_cmd = 'cp third-party/gfakluge/libgfakluge.a lib/'
    cp_hpp_cmd = 'mkdir -p include/gfakluge && cp third-party/gfakluge/src/*.hpp include/gfakluge'

    return {'title':    title_with_actions,
            'actions':  [build_cmd, cp_lib_cmd, cp_hpp_cmd],
            'targets':  ['lib/libgfakluge.a'],
            'uptodate': [run_once],
            'task_dep': ['display_libboink_config'],
            'clean':    [clean_targets,
                         (clean_folder, ['include/gfakluge'])]}


def task_deploy_cqf():
    build_cmd  = 'cd third-party/cqf && make'
    cp_hpp_cmd = 'mkdir -p include/cqf && cp third-party/cqf/gqf.h include/cqf/'
    cp_obj_cmd = 'mkdir -p {}/ && cp third-party/cqf/gqf.o {}/'.format(LIBBUILDDIR, LIBBUILDDIR)

    return {'title': title_with_actions,
            'actions': [build_cmd, cp_hpp_cmd, cp_obj_cmd],
            'targets': ['{}/gqf.o'.format(LIBBUILDDIR)],
            'file_dep': ['third-party/cqf/gqf.h'],
            'task_dep': ['display_libboink_config'],
            'clean': [clean_targets]}


def task_deploy_smhasher():
    build_cmd = [CXX, '-Ithird-party/smhasher/'] + CXXFLAGS + \
                ['-c', 'third-party/smhasher/MurmurHash3.cc', '-o', 'third-party/smhasher/MurmurHash3.o']
    cp_hpp_cmd = 'mkdir -p include/smhasher && cp third-party/smhasher/MurmurHash3.h include/smhasher/'
    cp_obj_cmd = 'mkdir -p {}/ && cp third-party/smhasher/MurmurHash3.o {}/'.format(LIBBUILDDIR, LIBBUILDDIR)

    return {'title': title_with_actions,
            'actions': [' '.join(build_cmd), cp_hpp_cmd, cp_obj_cmd],
            'targets': ['{}/MurmurHash3.o'.format(LIBBUILDDIR)],
            'file_dep': ['third-party/smhasher/MurmurHash3.h'],
            'task_dep': ['display_libboink_config'],
            'clean': [clean_targets]}


def task_deploy_sparsepp():

    return {'title':    title_with_actions,
            'actions':  ['cp -r third-party/sparsepp/sparsepp include/'],
            'task_dep': ['display_libboink_config'],
            'uptodate': [run_once],
            'clean':    [(clean_folder, ['include/sparsepp'])]}


def task_prepare_thirdparty():
    return {'actions':  None,
            'task_dep': ['display_libboink_config'],
            'task_dep': ['deploy_prometheus',
                         'deploy_sparsepp',
                         'deploy_smhasher',
                         'deploy_cqf',
                         'deploy_gfakluge',
                         'create_libboink_build_dirs']}

#
# Compile libboink object files
#
@create_after('prepare_thirdparty')
def task_compile_libboink():
    for source, obj in zip(SOURCE_FILES, OBJECT_FILES):
        source_file = os.path.join(*source)
        object_file = os.path.join(*obj)

        yield {'name':     object_file,
               'title':    title_with_actions,
               'file_dep': get_gcc_includes(source_file, INCLUDES, PKG) + [source_file],
               'task_dep': ['display_libboink_config'],
               'targets':  [object_file],
               'actions':  [cxx_command(source_file, object_file)],
               'clean':    True}


@create_after('compile_libboink')
def task_link_libboink():
    objects = [os.path.join(*obj) for obj in OBJECT_FILES] \
              + ['{}/gqf.o'.format(LIBBUILDDIR), '{}/MurmurHash3.o'.format(LIBBUILDDIR)]
    objects.sort(key=lambda path: os.path.basename(path))
    link_action = link_command(objects, LIBBOINKSO)
    ln_action   = ' '.join(['ln -sf',
                            os.path.abspath(LIBBOINKSO),
                            os.path.abspath(LIBBOINKSHARED)])

    return {'title':    title_with_actions,
            'file_dep': objects,
            'task_dep': ['display_libboink_config'],
            'targets':  [LIBBOINKSO,
                         LIBBOINKSHARED],
            'actions':  [link_action, ln_action],
            'clean':    True}


def task_boink_pc():
    src    = os.path.join(SOURCE_DIR, '{0}.pc.in'.format(PKG))
    target = os.path.join(LIBBUILDDIR, replace_ext(os.path.basename(src), ''))
    cmd = ("sed -e 's,@prefix@,{prefix},'  "
           "-e 's,@VERSION@,{version},' {src} >{dst}".format(prefix=PREFIX,
                                                             version=VERSION,
                                                             src=src,
                                                             dst=target))

    return {'title':    title_with_actions,
            'file_dep': [src],
            'task_dep': ['display_libboink_config'],
            'targets':  [target],
            'actions':  [cmd],
            'clean':    True}


@create_after('compile_libboink')
def task_boink_ranlib():
    target  = os.path.join(LIBBUILDDIR, 'lib{0}.a'.format(PKG))
    objects = [os.path.join(*obj) for obj in OBJECT_FILES]

    return {'title':    title_with_actions,
            'file_dep': objects,
            'task_dep': ['display_libboink_config'],
            'targets':  [target],
            'actions':  ['ar rcs {0} {1}'.format(target, ' '.join(objects)),
                         'ranlib {0}'.format(target)],
            'clean':    True}

def task_libboink():
    return {'actions': None,
            'task_dep': ['link_libboink',
                         'compile_libboink',
                         'boink_pc',
                         'boink_ranlib']}


def setupcmd(cmd):
    return ' '.join(['python', 'setup.py'] + cmd)


def task_install():
    return {'title': title_with_actions,
            'actions': [setupcmd(['install'])],
            'uptodate': [False]}


def task_test():
    return {'actions': [setupcmd(['develop']),
                        'py.test --pyargs boink.tests']}



