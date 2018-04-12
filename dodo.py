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

from jinja2 import Environment, PackageLoader, select_autoescape

from build_utils import (check_for_openmp,
                         distutils_dir_name,
                         build_dir,
                         BoinkReporter,
						 title_with_actions,
                         replace_ext,
                         resolve_dependencies)
from boink.types import PARAMS

DOIT_CONFIG  = {'verbosity': 2,
                'reporter': BoinkReporter,
                'default_tasks': ['build']}


PKG          = 'boink'
VERSION      = open(os.path.join(PKG, 'VERSION')).read().strip()
MOD_EXT      = sysconfig.get_config_var('SO')
PYX_FILES    = list(glob.glob(os.path.join(PKG, '*.pyx')))
PYX_NAMES    = [os.path.basename(fn) for fn in PYX_FILES]
PXD_FILES    = list(glob.glob(os.path.join(PKG, '*.pxd')))
PXD_NAMES    = [os.path.basename(fn) for fn in PXD_FILES]
PY_FILES     = list(glob.glob('boink/**/*.py', recursive=True))
EXT_SOS      = [pyx[:-4]+MOD_EXT for pyx in PYX_FILES]

INCLUDE_DIR  = os.path.join('include', PKG)
HEADERS      = ['dbg.hh', 'cdbg.hh', 'hashing.hh', 'assembly.hh', 'boink.hh',
                'consumer.hh', 'minimizers.hh', 'boink.hh']
HEADER_NAMES = HEADERS
SOURCES      = [filename[:-3] + '.cc' for filename in HEADERS]
SOURCE_NAMES = SOURCES
SOURCES      = [os.path.join('src', PKG, filename) for filename in SOURCES]
HEADERS      = [os.path.join('include', 'boink', filename) for filename in HEADERS]
SOURCES      = [filename for filename in SOURCES if os.path.isfile(filename)]
OBJECTS      = [os.path.splitext(source_file)[0] + '.o' for source_file in SOURCES]


DEP_MAP = resolve_dependencies(PYX_FILES, PYX_NAMES, PXD_FILES, PXD_NAMES,
                               SOURCES, SOURCE_NAMES, HEADERS, HEADER_NAMES,
                               EXT_SOS, MOD_EXT, INCLUDE_DIR, PKG)
pprint(DEP_MAP)
# Start libboink compile vars

PROFILING   = False
PROFILER    = 'gprof'
COLOR       = True
DEBUG       = False
DEBUG_FLAGS = ['-g']
PREFIX      = '/usr/local'

CXX         = os.environ.get('CXX', 'cc')
INCLUDES    = ['-I', os.path.abspath('include/')]
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
    CXXFLAGS += ['-arch',
                 'x86_64',
                 '-mmacosx-version-min=10.7',
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


def task_display_libboink_config():

    def print_config():
        print('CXX:', CXX)
        print('SOURCES:', *SOURCES)
        print('HEADERS:', *HEADERS)
        print('OBJECTS:', *OBJECTS)
        print('Compile:', cxx_command('SOURCE.cc', 'TARGET.o'))
        print('Link:', link_command(['OBJECT_1.o', 'OBJECT_2.o',
                                     'OBJECT_N.o'], SONAME))
    return {'actions': [print_config],
            'uptodate': [False]}


def task_compile_libboink():
    for source_file in SOURCES:
        object_file = os.path.splitext(source_file)[0] + '.o'

        yield {'name': object_file,
               'title': title_with_actions,
               'file_dep': HEADERS + [source_file],
               'task_dep': ['display_libboink_config'],
               'targets': [object_file],
               'actions': [cxx_command(source_file, object_file)],
               'clean': True}


def task_link_libboink():
    objects     = [os.path.splitext(fn)[0] + '.o' for fn in SOURCES]
    link_action = link_command(objects, LIBBOINKSO)
    ln_action   = ' '.join(['ln -sf',
                            os.path.abspath(LIBBOINKSO),
                            os.path.abspath(LIBBOINKSHARED)])

    return {'title': title_with_actions,
            'file_dep': objects,
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
            'file_dep': OBJECTS,
            'task_dep': ['display_libboink_config'],
            'targets': [target],
            'actions': ['ar rcs {0} {1}'.format(target, ' '.join(OBJECTS)),
                        'ranlib {0}'.format(target)],
            'clean': True}

def task_libboink():
    return {'actions': None,
            'task_dep': ['link_libboink',
                         'compile_libboink',
                         'boink_pc',
                         'boink_ranlib']}


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
        result['params'] = ','.join((p.strip('_') for p in param_bundle))
        results.append(result)
        cppclasses.append('{0}[{1}]'.format(cppclass, result['params']))
    return results, cppclasses


def get_templates(mod_prefix):
    env = get_template_env()
    
    pxd_template = env.get_template('{0}.tpl.pxd'.format(mod_prefix))
    pxd_dst_path = os.path.join('boink', '{0}.pxd.pxi'.format(mod_prefix))
    pyx_template = env.get_template('{0}.tpl.pyx'.format(mod_prefix))
    pyx_dst_path = os.path.join('boink', '{0}.pyx.pxi'.format(mod_prefix))

    return (pxd_template, pxd_dst_path,
            pyx_template, pyx_dst_path)


def render(pxd_template, pxd_dst_path,
           pyx_template, pyx_dst_path, type_bundles):
    
    with open(pxd_dst_path, 'w') as fp:
        rendered = pxd_template.render(dst_filename=pxd_dst_path,
                                       tpl_filename=pxd_template.name,
                                       type_bundles=type_bundles)
        fp.write(rendered)
    with open(pyx_dst_path, 'w') as fp:
        rendered = pyx_template.render(dst_filename=pyx_dst_path,
                                       tpl_filename=pyx_template.name,
                                       type_bundles=type_bundles)
        fp.write(rendered)


def setupcmd(cmd):
    return ' '.join(['python', 'setup.py'] + cmd)


def jinja_render_task(mod_prefix, type_dict):

    pxd_tpl, pxd_dst, pyx_tpl, pyx_dst = get_templates(mod_prefix)
    type_bundles, _ = generate_cpp_params(type_dict)
                                                                 

    return {'name': '{0}_types'.format(mod_prefix),
            'actions': [(render, [pxd_tpl, pxd_dst,
                                  pyx_tpl, pyx_dst,
                                  type_bundles])],
            'targets': [pxd_dst, pyx_dst],
            'file_dep': [os.path.join(PKG, 'templates', pxd_tpl.name),
                         os.path.join(PKG, 'templates', pyx_tpl.name)]}


def task_render_extension_classes():
    dbg_types = OrderedDict(storage_type=PARAMS['StorageTypes'],
                            shifter_type=PARAMS['ShifterTypes'])
    yield jinja_render_task('dbg', dbg_types)
    yield jinja_render_task('assembly', dbg_types)


def task_install():
    return {'actions': [setupcmd(['install'])]}


def task_test():
    return {'actions': [setupcmd(['develop']),
                        'py.test --pyargs boink.tests']}


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


@create_after(executed='render_extension_classes',
              target_regex='.*\\'+MOD_EXT)
def task_build():
    return {'file_dep': list(glob.glob('boink/*.pxi')) +
                        PYX_FILES +
                        PXD_FILES +
                        HEADERS +
                        SOURCES,
            'targets': EXT_SOS,
            'actions': [setupcmd(['build_ext', '--inplace'])],
            'clean':   CLEAN_ACTIONS}

