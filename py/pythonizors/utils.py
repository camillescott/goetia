import re


def match_template(short_name, full_name):
    '''
    Check that a given full name the template instantiation for the template named by short_name.
    '''
    expr = re.compile(r'{0}<[\S\s]*>'.format(short_name))
    return expr.match(full_name)


def is_template_inst(short_name, full_name):
    match = match_template(short_name.strip(), full_name.strip())
    if not match:
        return False
    else:
        return match.string == full_name.strip()
