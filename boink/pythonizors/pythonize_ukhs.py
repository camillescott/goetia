from numpy import zeros
import json

from boink.pythonizors import utils

def pythonize_boink_signatures(klass, name):
    is_inst, template =  utils.is_template_inst(name, 'UnikmerSignature')
    if is_inst:
        def signature(self):
            sig = zeros(len(self))
            raw = self.get_signature()
            for i, count in enumerate(raw):
                sig[i] = count
            return sig

        def __len__(self):
            return self.get_size()

        def to_dict(self, name):
            data = {'signature': self.signature.tolist(),
                    'W'        : self.K,
                    'K'        : self.bucket_K,
                    'size'     : len(self),
                    'name'     : name}
            return data

        def save(self, stream, name):
            data = [self.to_dict(name)]

            json.dump(data, stream)

        klass.Signature.signature = property(signature)
        klass.Signature.__len__   = __len__
        klass.Signature.to_dict   = to_dict
        klass.Signature.save      = save
