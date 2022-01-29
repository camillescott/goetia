from goetia import libgoetia

class cDBGSerialization:

    FORMATS = ['graphml', 'edgelist',
               'gfa1', 'gfa2', 'fasta', 'gml']

    def enum_from_str(file_format):
        if file_format in cDBGSerialization.FORMATS:
            if file_format == 'graphml':
                return libgoetia.GRAPHML
            elif file_format == 'fasta':
                return libgoetia.FASTA
            elif file_format == 'gfa1':
                return libgoetia.GFA1
            else:
                raise NotImplementedError("Support for {0} not yet "
                                          "implemented".format(file_format))
        else:
            formats = ', '.join(cDBGSerialization.FORMATS)
            raise ValueError("{0} not a valid save format. "
                             "Format must be one of: {1}".format(file_format,
                                                                 formats))
