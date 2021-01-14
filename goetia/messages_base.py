"""Module generated by SchemaModuleGenerator"""

from schemapi import SchemaBase, Undefined


class Messages(SchemaBase):
    """Messages schema wrapper

    oneOf(:class:`Interval`, :class:`SampleStarted`, :class:`SampleFinished`,
    :class:`SampleSaturated`, :class:`Error`, :class:`DistanceCalc`, :class:`EndStream`)
    """
    _schema = {'definitions': {'DistanceCalc': {'properties': {'distance': {'maximum': 1.0,
                                                                  'minimum': 0.0,
                                                                  'type': 'number'},
                                                     'file_names': {'type': 'array'},
                                                     'msg_type': {'enum': ['DistanceCalc'],
                                                                  'type': 'string'},
                                                     'sample_name': {'type': 'string'},
                                                     'sequence': {'minimum': 0,
                                                                  'type': 'integer'},
                                                     'stat': {'type': 'number'},
                                                     'stat_type': {'type': 'string'},
                                                     't': {'minimum': 0,
                                                           'type': 'integer'}},
                                      'required': ['msg_type',
                                                   't',
                                                   'sample_name',
                                                   'distance',
                                                   'stat',
                                                   'stat_type',
                                                   'file_names'],
                                      'type': 'object'},
                     'EndStream': {'properties': {'msg_type': {'enum': ['EndStream'],
                                                               'type': 'string'}},
                                   'required': ['msg_type'],
                                   'type': 'object'},
                     'Error': {'properties': {'error': {'type': 'string'},
                                              'file_names': {'type': 'array'},
                                              'msg_type': {'enum': ['Error'],
                                                           'type': 'string'},
                                              'sample_name': {'type': 'string'},
                                              'sequence': {'minimum': 0,
                                                           'type': 'integer'},
                                              't': {'minimum': 0,
                                                    'type': 'integer'}},
                               'required': ['msg_type',
                                            't',
                                            'sample_name',
                                            'error',
                                            'file_names'],
                               'type': 'object'},
                     'Interval': {'properties': {'file_names': {'type': 'array'},
                                                 'modulus': {'minimum': 0,
                                                             'type': 'integer'},
                                                 'msg_type': {'const': 'Interval',
                                                              'default': 'Interval',
                                                              'type': 'string'},
                                                 'sequence': {'minimum': 0,
                                                              'type': 'integer'},
                                                 't': {'minimum': 0,
                                                       'type': 'integer'}},
                                  'required': ['msg_type',
                                               't',
                                               'sample_name',
                                               'file_names'],
                                  'type': 'object'},
                     'SampleFinished': {'properties': {'file_names': {'type': 'array'},
                                                       'msg_type': {'enum': ['SampleFinished'],
                                                                    'type': 'string'},
                                                       'sample_name': {'type': 'string'},
                                                       'sequence': {'minimum': 0,
                                                                    'type': 'integer'},
                                                       't': {'minimum': 0,
                                                             'type': 'integer'}},
                                        'required': ['msg_type',
                                                     'sample_name',
                                                     't',
                                                     'file_names'],
                                        'type': 'object'},
                     'SampleSaturated': {'properties': {'file_names': {'type': 'array'},
                                                        'msg_type': {'enum': ['SampleSaturated'],
                                                                     'type': 'string'},
                                                        'sample_name': {'type': 'string'},
                                                        'sequence': {'minimum': 0,
                                                                     'type': 'integer'},
                                                        't': {'minimum': 0,
                                                              'type': 'integer'}},
                                         'required': ['msg_type',
                                                      'sample_name',
                                                      't',
                                                      'file_names'],
                                         'type': 'object'},
                     'SampleStarted': {'properties': {'file_names': {'type': 'array'},
                                                      'msg_type': {'enum': ['SampleStarted'],
                                                                   'type': 'string'},
                                                      'sample_name': {'type': 'string'}},
                                       'required': ['msg_type',
                                                    'sample_name',
                                                    'file_names'],
                                       'type': 'object'}},
     'oneOf': [{'$ref': '#/definitions/Interval'},
               {'$ref': '#/definitions/SampleStarted'},
               {'$ref': '#/definitions/SampleFinished'},
               {'$ref': '#/definitions/SampleSaturated'},
               {'$ref': '#/definitions/Error'},
               {'$ref': '#/definitions/DistanceCalc'},
               {'$ref': '#/definitions/EndStream'}]}
    _rootschema = _schema

    def __init__(self, *args, **kwds):
        super(Messages, self).__init__(*args, **kwds)



class Interval(SchemaBase):
    """Interval schema wrapper

    Mapping(required=[msg_type, t, sample_name, file_names])

    Attributes
    ----------

    file_names : List(Mapping(required=[]))

    msg_type : string

    t : integer

    modulus : integer

    sequence : integer

    """
    _schema = {'$ref': '#/definitions/Interval'}
    _rootschema = Messages._schema

    def __init__(self, file_names=Undefined, msg_type=Undefined, sample_name=Undefined, t=Undefined,
                 modulus=Undefined, sequence=Undefined, **kwds):
        super(Interval, self).__init__(file_names=file_names, msg_type=msg_type,
                                       sample_name=sample_name, t=t, modulus=modulus, sequence=sequence,
                                       **kwds)



class SampleStarted(SchemaBase):
    """SampleStarted schema wrapper

    Mapping(required=[msg_type, sample_name, file_names])

    Attributes
    ----------

    file_names : List(Mapping(required=[]))

    msg_type : enum('SampleStarted')

    sample_name : string

    """
    _schema = {'$ref': '#/definitions/SampleStarted'}
    _rootschema = Messages._schema

    def __init__(self, file_names=Undefined, msg_type=Undefined, sample_name=Undefined, **kwds):
        super(SampleStarted, self).__init__(file_names=file_names, msg_type=msg_type,
                                            sample_name=sample_name, **kwds)



class SampleFinished(SchemaBase):
    """SampleFinished schema wrapper

    Mapping(required=[msg_type, sample_name, t, file_names])

    Attributes
    ----------

    file_names : List(Mapping(required=[]))

    msg_type : enum('SampleFinished')

    sample_name : string

    t : integer

    sequence : integer

    """
    _schema = {'$ref': '#/definitions/SampleFinished'}
    _rootschema = Messages._schema

    def __init__(self, file_names=Undefined, msg_type=Undefined, sample_name=Undefined, t=Undefined,
                 sequence=Undefined, **kwds):
        super(SampleFinished, self).__init__(file_names=file_names, msg_type=msg_type,
                                             sample_name=sample_name, t=t, sequence=sequence, **kwds)



class SampleSaturated(SchemaBase):
    """SampleSaturated schema wrapper

    Mapping(required=[msg_type, sample_name, t, file_names])

    Attributes
    ----------

    file_names : List(Mapping(required=[]))

    msg_type : enum('SampleSaturated')

    sample_name : string

    t : integer

    sequence : integer

    """
    _schema = {'$ref': '#/definitions/SampleSaturated'}
    _rootschema = Messages._schema

    def __init__(self, file_names=Undefined, msg_type=Undefined, sample_name=Undefined, t=Undefined,
                 sequence=Undefined, **kwds):
        super(SampleSaturated, self).__init__(file_names=file_names, msg_type=msg_type,
                                              sample_name=sample_name, t=t, sequence=sequence, **kwds)



class Error(SchemaBase):
    """Error schema wrapper

    Mapping(required=[msg_type, t, sample_name, error, file_names])

    Attributes
    ----------

    error : string

    file_names : List(Mapping(required=[]))

    msg_type : enum('Error')

    sample_name : string

    t : integer

    sequence : integer

    """
    _schema = {'$ref': '#/definitions/Error'}
    _rootschema = Messages._schema

    def __init__(self, error=Undefined, file_names=Undefined, msg_type=Undefined, sample_name=Undefined,
                 t=Undefined, sequence=Undefined, **kwds):
        super(Error, self).__init__(error=error, file_names=file_names, msg_type=msg_type,
                                    sample_name=sample_name, t=t, sequence=sequence, **kwds)



class DistanceCalc(SchemaBase):
    """DistanceCalc schema wrapper

    Mapping(required=[msg_type, t, sample_name, distance, stat, stat_type, file_names])

    Attributes
    ----------

    distance : float

    file_names : List(Mapping(required=[]))

    msg_type : enum('DistanceCalc')

    sample_name : string

    stat : float

    stat_type : string

    t : integer

    sequence : integer

    """
    _schema = {'$ref': '#/definitions/DistanceCalc'}
    _rootschema = Messages._schema

    def __init__(self, distance=Undefined, file_names=Undefined, msg_type=Undefined,
                 sample_name=Undefined, stat=Undefined, stat_type=Undefined, t=Undefined,
                 sequence=Undefined, **kwds):
        super(DistanceCalc, self).__init__(distance=distance, file_names=file_names, msg_type=msg_type,
                                           sample_name=sample_name, stat=stat, stat_type=stat_type, t=t,
                                           sequence=sequence, **kwds)



class EndStream(SchemaBase):
    """EndStream schema wrapper

    Mapping(required=[msg_type])

    Attributes
    ----------

    msg_type : enum('EndStream')

    """
    _schema = {'$ref': '#/definitions/EndStream'}
    _rootschema = Messages._schema

    def __init__(self, msg_type=Undefined, **kwds):
        super(EndStream, self).__init__(msg_type=msg_type, **kwds)
