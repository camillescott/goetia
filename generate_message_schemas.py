import schemapi

msg_schema = {
    'definitions': {
        'Interval': {
            "type": "object",
            "required": ["msg_type", "t", "state", 'sample_name'],
            'properties': {
                'msg_type': {'type': 'string',
                             'const': 'Interval',
                             'default': 'Interval'},
                't': {'type': 'integer',
                       'minimum': 0},
                'state': {'type': 'string',
                          'enum': ['fine', 'medium', 'coarse']}
            }
        },
        'SampleStarted': {
            "type": "object",
            "required": ["msg_type", "sample_name"],
            "properties": {
                "msg_type": {"type": "string",
                             "enum": ['SampleStarted']},
                "sample_name": {"type": "string"}
            }
        },
        "SampleFinished": {
            "type": "object",
            "required": ["msg_type", "sample_name", "t"],
            "properties": {
                'msg_type': {'type': 'string',
                             'enum': ['SampleFinished']},
                't': {'type': 'integer',
                      'minimum': 0},
                "sample_name": {"type": "string"}
            }
        },
        "Error": {
            "type": 'object',
            'required': ['msg_type', 't', 'sample_name', 'error'],
            'properties': {
                'msg_type': {'type': 'string',
                              'enum': ['Error']},
                't': {'type': 'integer',
                      'minimum': 0},
                "sample_name": {"type": "string"},
                "error": {"type": "string"}
            }
        },
        "DistanceCalc": {
            "type": 'object',
            'required': ['msg_type', 't', 'sample_name', 'delta', 'distance'],
            'properties': {
                'msg_type': {'type': 'string',
                             'enum': ['DistanceCalc']},
                't': {'type': 'integer',
                      'minimum': 0},
                'sample_name': {'type': 'string'},
                'delta': {'type': 'integer',
                          'minimum': 0},
                'distance': {'type': 'number',
                             'minimum': 0.0,
                             'maximum': 1.0}
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
    "oneOf": [
        {"$ref": "#/definitions/Interval"},
        {"$ref": "#/definitions/SampleStarted"},
        {"$ref": "#/definitions/SampleFinished"},
        {"$ref": "#/definitions/Error"},
        {"$ref": "#/definitions/DistanceCalc"},
        {"$ref": "#/definitions/EndStream"}
    ]
}

api = schemapi.SchemaModuleGenerator(msg_schema, root_name='Messages')
api.write_module('boink/messages_base.py')
