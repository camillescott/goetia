from goetia import schemapi

msg_schema = {
    'definitions': {
        'Interval': {
            'type': 'object',
            'required': ['msg_type', 't', 'state', 'sample_name', 'file_names'],
            'properties': {
                'msg_type': {'type': 'string',
                             'const': 'Interval',
                             'default': 'Interval'},
                't': {'type': 'integer',
                       'minimum': 0},
                'state': {'type': 'string',
                          'enum': ['fine', 'medium', 'coarse']},
                'file_names': {'type': 'array'}
            }
        },
        'SampleStarted': {
            'type': 'object',
            'required': ['msg_type', 'sample_name', 'file_names'],
            'properties': {
                'msg_type': {'type': 'string',
                             'enum': ['SampleStarted']},
                'sample_name': {'type': 'string'},
                'file_names': {'type': 'array'}
            }
        },
        'SampleFinished': {
            'type': 'object',
            'required': ['msg_type', 'sample_name', 't', 'file_names'],
            'properties': {
                'msg_type': {'type': 'string',
                             'enum': ['SampleFinished']},
                't': {'type': 'integer',
                      'minimum': 0},
                'sample_name': {'type': 'string'},
                'file_names': {'type': 'array'}
            }
        },
        'SampleSaturated': {
            'type': 'object',
            'required': ['msg_type', 'sample_name', 't', 'file_names'],
            'properties': {
                'msg_type': {'type': 'string',
                             'enum': ['SampleSaturated']},
                't': {'type': 'integer',
                      'minimum': 0},
                'sample_name': {'type': 'string'},
                'file_names': {'type': 'array'}
            }
        },
        'Error': {
            'type': 'object',
            'required': ['msg_type', 't', 'sample_name', 'error', 'file_names'],
            'properties': {
                'msg_type': {'type': 'string',
                              'enum': ['Error']},
                't': {'type': 'integer',
                      'minimum': 0},
                'sample_name': {'type': 'string'},
                'error': {'type': 'string'},
                'file_names': {'type': 'array'}
            }
        },
        'DistanceCalc': {
            'type': 'object',
            'required': ['msg_type', 't', 'sample_name', 'delta', 'distance', 'stat', 'stat_type', 'file_names'],
            'properties': {
                'msg_type': {'type': 'string',
                             'enum': ['DistanceCalc']},
                't': {'type': 'integer',
                      'minimum': 0},
                'sample_name': {'type': 'string'},
                'stat_type': {'type': 'string'},
                'delta': {'type': 'integer',
                          'minimum': 0},
                'distance': {'type': 'number',
                             'minimum': 0.0,
                             'maximum': 1.0},
                'stat': {'type': 'number'},
                'file_names': {'type': 'array'}
            }
        },
        'EndStream': {
            'type': 'object',
            'required': ['msg_type'],
            'properties': {
                'msg_type': {'type': 'string',
                              'enum': ['EndStream']}
            }
        }

    },
    'oneOf': [
        {'$ref': '#/definitions/Interval'},
        {'$ref': '#/definitions/SampleStarted'},
        {'$ref': '#/definitions/SampleFinished'},
        {'$ref': '#/definitions/SampleSaturated'},
        {'$ref': '#/definitions/Error'},
        {'$ref': '#/definitions/DistanceCalc'},
        {'$ref': '#/definitions/EndStream'}
    ]
}

api = schemapi.SchemaModuleGenerator(msg_schema, root_name='Messages')
api.write_module('goetia/messages_base.py')
