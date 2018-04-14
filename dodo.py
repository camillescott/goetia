#!/usr/bin/env doit

from collections import OrderedDict
import glob
import itertools
import os
from pprint import pprint
import sysconfig
import sys

from doit import create_after
from doit.task import clean_targets
import yaml

from build_utils import (check_for_openmp,
                         distutils_dir_name,
                         build_dir,
                         lib_dir,
                         BoinkReporter,
						 title_with_actions,
                         replace_ext,
                         replace_exts,
                         resolve_dependencies,
                         isfile_filter,
                         flatten,
                         render,
                         get_templates,
                         generate_cpp_params)

DOIT_CONFIG  = {'verbosity': 2,
                'reporter': BoinkReporter,
                'default_tasks': ['extensions']}
EXT_META     = yaml.load(open('extensions.yaml').read())


PKG          = 'boink'
VERSION      = open(os.path.join(PKG, 'VERSION')).read().strip()
MOD_EXT      = sysconfig.get_config_var('SO')
MOD_NAMES    = [m['name'] for m in EXT_META['extensions']]
PYX_FILES    = {m: os.path.join(PKG, '{0}.pyx'.format(m)) for m in MOD_NAMES}
PXD_FILES    = isfile_filter(replace_exts(PYX_FILES, '.pxd'))
TPL_FILES    = {m['name']:m['templates'] for m in EXT_META['extensions'] \
                if 'templates' in m} 
PY_FILES     = list(glob.glob('boink/**/*.py', recursive=True))
MOD_FILES    = {m:m+MOD_EXT for m in MOD_NAMES}

INCLUDE_DIR  = os.path.join('include', PKG)
SOURCE_DIR   = os.path.join('src', PKG)
HEADER_NAMES = EXT_META['headers']
HEADER_FILES = [os.path.join(INCLUDE_DIR, header) for header in HEADER_NAMES]
SOURCE_NAMES = replace_exts(HEADER_NAMES, '.cc')
SOURCE_FILES = isfile_filter([os.path.join(SOURCE_DIR, source) \
                             for source in SOURCE_NAMES])
OBJECT_FILES = replace_exts(SOURCE_FILES, '.o')

DEP_MAP = resolve_dependencies(PYX_FILES,
                               PXD_FILES,
                               MOD_FILES,
                               MOD_EXT,
                               INCLUDE_DIR,
                               PKG)
#pprint(PYX_FILES)
#pprint(PXD_FILES)
#pprint(TPL_FILES)
#pprint(DEP_MAP)


# Start libboink compile vars

PROFILING   = False
PROFILER    = 'gprof'
COLOR       = True
DEBUG       = False
DEBUG_FLAGS = ['-g']
PREFIX      = '/usr/local'

CXX         = os.environ.get('CXX', 'cc')
INCLUDES    = ['-I', os.path.abspath('include/'), '-I.']
WARNINGS    = ['-Wall']
COMMON      = ['-O3', '-fPIC', '-fno-omit-frame-pointer']

CPPFLAGS    = []

CXXFLAGS    = COMMON + WARNINGS
CXXFLAGS   += ['-Wstrict-null-sentinel', '-std=c++14']
CXXFLAGS   += INCLUDES
CXXFLAGS   += CPPFLAGS

CFLAGS      = ['-Wshadow', '-Wcast-align', '-Wstrict-prototypes']
CFLAGS     += INCLUDES
CFLAGS     += CPPFLAGS

LDFLAGS     = ['-loxli']

if DEBUG:
    CXXFLAGS += DEBUG_FLAGS
    CFLAGS   += DEBUG_FLAGS

if PROFILER == 'gprof':
    CXXFLAGS += ['-pg']
    CFLAGS   += ['-pg']
    LDFLAGS  += ['-pg']

if sys.platform == 'linux':
    LDFLAGS  += ['-pthread']

if check_for_openmp():
    LDFLAGS  += ['-fopenmp']
    CXXFLAGS += ['-fopenmp']

if COLOR:
    CXXFLAGS += ['-fdiagnostics-color']

if sys.platform == 'darwin':
    CXXFLAGS += [ '-mmacosx-version-min=10.7',
                 '-stdlib=libc++']

LIB_VERSION = '1'
LIBDIR = os.path.join('src', PKG)

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

LIBBOINKSO     = os.path.join(LIBDIR, SONAME)
LIBBOINKSHARED = os.path.join(LIBDIR, 'lib{0}.{1}'.format(PKG, SHARED_EXT))
CXXFLAGS  += ['-DVERSION={0}'.format(VERSION)]


CYTHON_DIRECTIVES = \
{
    'embedsignature': True,
    'c_string_type': 'unicode',
    'c_string_encoding': 'utf8'
}
CYTHON_FLAGS     = ['-X', ','.join(('{0}={1}'.format(k, v) \
                                    for (k, v) in CYTHON_DIRECTIVES.items())),
                    '-3',
                    '--line-directives',
                    '--cplus']

new_types = []
for type_dict in EXT_META['types']:
    if type_dict['types'] is not None:
        type_dict['types'] = ['_'+_type for _type in type_dict['types']]
    #new_types.append(type_dict)
#EXT_META['types'] = new_types
pprint(EXT_META)


def cxx_command(source, target):
    cmd = [ CXX ] \
          + CXXFLAGS \
          + LDFLAGS \
          + ['-c -o', target, source]
    return ' '.join(cmd)


def link_command(objects, target):
    cmd = [ CXX ] \
          + CXXFLAGS \
          + LDFLAGS \
          + ['-shared', '-o', target] \
          + objects
    return ' '.join(cmd)


def cython_command(pyx_file):
    cmd = [ 'cython' ] \
          + ['-I include/'] \
          + CYTHON_FLAGS \
          + [pyx_file]
    return ' '.join(cmd)


def cy_cxx_command(cy_source, cy_dst):
    cmd = [ CXX ] \
          + [sysconfig.get_config_var('CFLAGS')] \
          + INCLUDES \
          + ['-I'+sysconfig.get_config_var('INCLUDEPY')] \
          + LDFLAGS \
          + ['-c ', cy_source] \
          + ['-o', cy_dst] \
          + CXXFLAGS
    return ' '.join(cmd)


def cy_link_command(cy_object, mod):
    PYLDSHARE = sysconfig.get_config_var('LDSHARED')
    PYLDSHARE = PYLDSHARE.split()[1:] # remove compiler invoke
    cmd = [ CXX ] \
          + PYLDSHARE \
          + [cy_object] \
          + LDFLAGS \
          + ['-o', mod]
    return ' '.join(cmd)


def task_display_libboink_config():

    def print_config():
        print('CXX:', CXX)
        print('SOURCES:', *SOURCE_FILES)
        print('HEADERS:', *HEADER_FILES)
        print('OBJECTS:', *OBJECT_FILES)
        print('Compile:', cxx_command('SOURCE.cc', 'TARGET.o'))
        print('Link:', link_command(['OBJECT_1.o', 'OBJECT_2.o',
                                     'OBJECT_N.o'], SONAME))
    return {'actions': [print_config],
            'uptodate': [False]}


def task_compile_libboink():
    for source_file, object_file in zip(SOURCE_FILES, OBJECT_FILES):

        yield {'name': object_file,
               'title': title_with_actions,
               'file_dep': HEADER_FILES + [source_file],
               'task_dep': ['display_libboink_config'],
               'targets': [object_file],
               'actions': [cxx_command(source_file, object_file)],
               'clean': True}


def task_link_libboink():
    link_action = link_command(OBJECT_FILES, LIBBOINKSO)
    ln_action   = ' '.join(['ln -sf',
                            os.path.abspath(LIBBOINKSO),
                            os.path.abspath(LIBBOINKSHARED)])

    return {'title': title_with_actions,
            'file_dep': OBJECT_FILES,
            'task_dep': ['display_libboink_config'],
            'targets': [LIBBOINKSO],
            'actions': [link_action, ln_action],
            'clean': True}


def task_boink_pc():
    target = os.path.join(LIBDIR, '{0}.pc'.format(PKG))
    src    = target + '.in'
    cmd = "sed -e 's,@prefix@,{prefix},'  "\
          "-e 's,@VERSION@,{version},' {src} >{dst}".format(prefix=PREFIX,
                                                            version=VERSION,
                                                            src=src,
                                                            dst=target)

    return {'title': title_with_actions,
            'file_dep': [src],
            'task_dep': ['display_libboink_config'],
            'targets': [target],
            'actions': [cmd],
            'clean': True}


def task_boink_ranlib():
    target = 'lib{0}.a'.format(PKG)

    return {'title': title_with_actions,
            'file_dep': OBJECT_FILES,
            'task_dep': ['display_libboink_config'],
            'targets': [target],
            'actions': ['ar rcs {0} {1}'.format(target, ' '.join(OBJECT_FILES)),
                        'ranlib {0}'.format(target)],
            'clean': True}

def task_libboink():
    return {'actions': None,
            'task_dep': ['link_libboink',
                         'compile_libboink',
                         'boink_pc',
                         'boink_ranlib']}

def setupcmd(cmd):
    return ' '.join(['python', 'setup.py'] + cmd)


def task_install():
    return {'actions': [setupcmd(['install'])]}


def task_test():
    return {'actions': [setupcmd(['develop']),
                        'py.test --pyargs boink.tests']}


def jinja_render_task(mod_prefix, type_dict):

    pxd_tpl, pxd_dst, pyx_tpl, pyx_dst = get_templates(mod_prefix)
    type_bundles, _ = generate_cpp_params(type_dict)
    print(type_bundles)

    return {'name': '{0}_types'.format(mod_prefix),
            'actions': [(render, [pxd_tpl, pxd_dst,
                                  pyx_tpl, pyx_dst,
                                  type_bundles])],
            'targets': [pxd_dst, pyx_dst],
            'file_dep': [os.path.join(PKG, 'templates', pxd_tpl.name),
                         os.path.join(PKG, 'templates', pyx_tpl.name)],
            'clean': True}


def task_render_extension_classes():
    types = {t['name']: t['types'] for t in EXT_META['types']}
    for mod_name, tpl_file in TPL_FILES.items():
        dbg_types = OrderedDict(storage_type=types['StorageType'],
                                shifter_type=types['ShifterType'])
        yield jinja_render_task(mod_name, dbg_types)


CLEAN_ACTIONS = \
[
    'rm -f {0}/*.cpp'.format(PKG),
    'rm -f {0}/*.so'.format(PKG),
    'find ./ -type d -name __pycache__ -exec rm -rf {} +',
    'find ./{0}/ -type f -name *{1} -exec rm -f {{}} +'.format(PKG, MOD_EXT),
    'find ./{0}/ -type f -name *.pyc -exec rm -f {{}} +'.format(PKG),
    'find ./{0}/ -type f -name *.cpp -exec rm -f {{}} +'.format(PKG),
    'rm -rf build dist {0}.egg-info'.format(PKG)
]


def task_cythonize():
    templated = flatten(list(TPL_FILES.values()))
    file_dep = [os.path.join(PKG, t)+'.pxi' for t in templated]
    for mod, source in PYX_FILES.items():
        target = replace_ext(source, '.cpp')
        yield { 'name': target,
                'title': title_with_actions,
                'file_dep': file_dep + [source],
                'targets': [target],
                'actions': ['mkdir -p ' + build_dir(),
                            cython_command(source)],
                'clean': True}


def task_compile_cython_cpp():
    for mod in MOD_NAMES:
        source = os.path.join(PKG, '{0}.cpp'.format(mod))
        target = os.path.join(build_dir(), '{0}.o'.format(mod))
        yield { 'name': target,
                'title': title_with_actions,
                'file_dep': [source],
                'targets': [target],
                'actions': [cy_cxx_command(source,
                                           target)],
                'clean': True}

def task_extensions():
    for mod, pyx_file in PYX_FILES.items():
        source = os.path.join(build_dir(),
                              '{0}.o'.format(mod))
        target = os.path.join(lib_dir(),
                              '{0}.cpython-36m-darwin.so'.format(mod))
        cp_target = os.path.join(PKG, os.path.basename(target))
        yield {'name': target,
                'title': title_with_actions,
                'file_dep': [source],
                'targets': [target, cp_target],
                'actions': ['mkdir -p ' + lib_dir(),
                            cy_link_command(source, target),
                            'cp {0} {1}'.format(target, cp_target)],
                'clean': True}
