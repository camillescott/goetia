import cppyy
import re


def match_template(full_name, short_name):
    '''
    Check that a given full name the template instantiation for the template named by short_name.
    '''
    expr = re.compile(r'{0}<[\S\s]*>'.format(short_name))
    match = expr.match(full_name)
    if match:
        template = full_name[slice(*match.span())]
        return match, template
    return match, ''


def is_template_inst(full_name, short_name):
    full_name = full_name.strip()
    match, template = match_template(full_name, short_name.strip())
    if not match:
        return False, None
    else:
        return full_name.startswith(template), template


def is_member(q_name, q_namespace):
    namespace, _, name = q_name.rpartition('::')
    return namespace == q_namespace, name


def convert_nullptr(func):
    
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        if cppyy.addressof(result) == 0:
            return None
        return result

    return wrapper
