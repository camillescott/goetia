import cppyy
from cppyy.gbl import std
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
    '''
    True if q_name is a member of q_namespace.
    '''
    namespace, _, name = q_name.rpartition('::')
    return namespace == q_namespace, name


def convert_nullptr(func):
    '''
    Decorator for functions that could return nullptr to return None.
    '''
    
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        if cppyy.addressof(result) == 0:
            return None
        return result

    return wrapper


def make_tuple_from_type(tt, *args):
    ''' SILLIEST FUNCTION EVER.

    Firstly, I can't figure out how to extract the template types directly from the cppyy template
    proxy for a typedef'd std::tuple. Secondly, you *have* to use std::make_tuple from Python-land
    to make properly-typed tuples that will type-match C++ function signatures expecting them. Thirdly,
    my variadic utility function defined at the c++ level which is meant to create a new tuple of
    the same type given an existing tuple doesn't work in Python-land without explicit declarations.
    Next, if I use `type` on a Python-land tuple type, it gets cast to a Python equivalentm so I have to
    use std::tuple_element. Finally, if I save the types in a list (as strings or otherwise) and then 
    pass them to the template proxy THEY GET EXTRACTED AND TURNED BACK INTO PYTHON TYPES, and I can't splat
    into the square brackets because it's invalid syntax. So, I'm left with... this. Goddamn, man.
    '''

    def name(tup_t, idx):
        return std.tuple_element[idx, tup_t].type.__cpp_name__

    if tt._tuple_len == 1:
        return std.make_tuple[name(tt, 0)](*args)
    elif tt._tuple_len == 2:
        return std.make_tuple[name(tt, 0), name(tt, 1)](*args)
    elif tt._tuple_len == 3:
        return std.make_tuple[name(tt, 0), name(tt, 1), name(tt, 2)](*args)
    elif tt._tuple_len == 4:
        return std.make_tuple[name(tt, 0), name(tt, 1), name(tt, 2), name(tt, 3)](*args)
    else:
        raise TypeError(f"Too many values in provided tuple (got {tt._tuple_len}, must be four or less).")