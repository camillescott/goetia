from goetia.pythonizors.utils import is_template_inst

def pythonize_goetia_parsing(klass, name):
    is_fastx, _ = is_template_inst(name, 'FastxParser')
    if is_fastx:
        def __iter__(self):
            while not self.is_complete():
                record = self.next()
                if record:
                    yield record
        klass.__iter__ = __iter__



    is_split, _ = is_template_inst(name, 'SplitPairedReader')
    if is_split:
        def __iter__(self):
            while not self.is_complete():
                pair = self.next()
                left = pair.left if pair.has_left else None
                right = pair.right if pair.has_right else None 
                if left is not None or right is not None:
                    yield left, right
        klass.__iter__ = __iter__