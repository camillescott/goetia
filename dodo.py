#!/usr/bin/env doit

from collections import OrderedDict
import glob
import itertools
import os

from doit import create_after
from doit.task import clean_targets
from jinja2 import Environment, PackageLoader, select_autoescape

from boink.types import PARAMS


PKG       = 'boink'
MODEXT    = sysconfig.get_config_var('SO')
PYX_FILES = list(glob.glob(os.path.join(PKG, '*.pyx')))
PXD_FILES = list(glob.glob(os.path.join(PKG, '*.pxd')))
PY_FILES  = list(glob.glob('boink/**/*.py', recursive=True))
EXT_SOS   = [pyx[:-4]+MOD_EXT for pyx in PYX_FILES]

HEADERS   = ['dbg.hh', 'cdbg.hh', 'hashing.hh', 'assembly.hh', 'boink.hh',
             'consumer.hh', 'minimizers.hh']
HEADERS   = [os.path.join('include', 'boink', filename) for filename in DEPENDS]

SOURCES   = [filename[:-2] + '.cc' for filename in DEPENDS]
SOURCES   = [filename for filename in SOURCES if os.path.isfile(filename)]


def get_template_env():
    return Environment(loader=PackageLoader('boink', 'templates'),
                       trim_blocks=True,
                       lstrip_blocks=True)


def generate_cpp_params(type_dict, cppclass='class'):
    results = []
    cppclasses = []
    for param_bundle in itertools.product(*(v for _, v in type_dict.items())):
        result = {}
        for name, param in zip(kwargs.keys(), param_bundle):
            result[name] = param
        result['suffix'] = '_'.join(param_bundle)
        result['params'] = ','.join((p.strip('_') for p in param_bundle))
        results.append(result)
        cppclasses.append('{0}[{1}]'.format(result['params']))
    return results, cppclasses


def get_templates(mod_prefix):
    env = get_template_env()
    
    pxd_template = env.get_template('{0}.tpl.pxd'.format(prefix))
    pxd_dst_path = os.path.join('boink', '{0}.pxd.pxi'.format(prefix))
    pyx_template = env.get_template('{0}.tpl.pyx'.format(prefix))
    pyx_dst_path = os.path.join('boink', '{0}.pyx.pxi'.format(prefix))

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
        res = pyx_tpl.render(dst_filename=pyx_dst_path,
                             tpl_filename=pyx_template.name,
                             type_list=type_list)
        fp.write(res)


def setupcmd(cmd):
    return ' '.join(['python', 'setup.py'] + cmd)


def jinja_render_task(mod_prefix, type_dict):

    pxd_tpl, pxd_dst, pyx_tpl, pyx_dst = get_template(mod_prefix)
    type_bundles, _ = generate_cpp_params(type_dict)
                                                                 

    return {'name': '{0}_types'.format(mod_prefix),
            'actions': [(render, [pxd_tpl, pxd_dst_path,
                                  pyx_tpl, pyx_dst_path,
                                  type_bundles])],
            'targets': [pxd_dst_path, pyx_dst_path],
            'file_dep': [pxd_tpl.name, pyx_tpl.name]}


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
            'clean': True}


def task_clean():
    actions = [
    'rm -f {0}/*.cpp'.format(PKG),
	'rm -f {0}/*.so'.format(PKG),
	'@find ./ -type d -name __pycache__ -exec rm -rf {} +',
	'@find ./{0}/ -type f -name *{1} -exec rm -f {{}} +'.format(PKG, MODEXT),
	'@find ./{0}/ -type f -name *.pyc -exec rm -f {{}} +'.format(PKG),
	'@find ./{0}/ -type f -name *.cpp -exec rm -f {{}} +'.format(PKG),
	'rm -rf build dist {0}.egg-info'.format(PKG)
    ]
    return {'actions': actions}

