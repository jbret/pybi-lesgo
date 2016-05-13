#!/usr/bin/python

import os
import struct
import numpy as np

def readmyfile(filename, bytes=8, endian='>d'):
    fileSize = os.path.getsize(filename)
    print "File size     : ", fileSize
    #print "3*nx*ny*nz2*8 : ", 3*nx*ny*nz2*8
    values = np.empty(fileSize/bytes)
    with open(filename, 'rb') as f:
        for i in range(len(values)):
            values[i] = struct.unpack(endian, f.read(bytes))[0]
	    #print "Val: ", i, values[i]
    return values

