from boink import boink as libboink


storage_types = [(libboink.storage.SparseppSetStorage, ()),
                 (libboink.storage.BitStorage, (100000, 4)),
                 (libboink.storage.ByteStorage, (100000, 4)),
                 (libboink.storage.NibbleStorage, (100000, 4))]

