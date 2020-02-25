from boink import messages_base

'''
Make the constructor pass in msg_type matching the
schema name as default. We dynamically create a subclass
with the new constructor and inject it into the globals()
for this module.

Eldritch horrors abound!
'''
def create_msg_type(base):
    def __init__(self, **kwargs):
        base.__init__(self, msg_type=base.__name__, **kwargs)
    globals()[base.__name__] = type(base.__name__, (base,), {'__init__': __init__})


'''
Now, iterator through all the members of `messages_base`
(which was generated by schemapi) and create the new
subclass for any attribute that is derived from
`SchemaBase` (ie, is a generated object), but that
isn't `SchemaBase` itself or the catch-all `Messages`
API.
'''
for attr in dir(messages_base):
    if attr in ('Messages', 'SchemaBase'):
        continue
    attr = getattr(messages_base, attr)
    if type(attr) is type and issubclass(attr, messages_base.SchemaBase):
        create_msg_type(attr)


'''
Dummy class that we'll use to represent all messages.
'''
class AllMessages:
    pass