"""
Reads lesgo's binary data output into a numpy array
Author: Joel Bretheim
"""
import numpy as np
from read_lesgo_bin import readmyfile

# neglect top point in each proc except last proc
# remember, 0:8 does not grab the point at index 8, just indices 0-7
def assemble_field(fileName, nproc, num, nz2, nz, ny, nx):
    nz_ = nz2 - 1
    out = np.zeros((num, nz, ny, nx))
    for i in range(0, nproc):
        filecontents = readmyfile( fileName + str(i) )
        out_i = np.reshape(filecontents, (num, nz2, ny, nx))
        a = i*nz_;  b = (i+1)*nz_
        out[:,a:b,:,:] = out_i[:,0:nz_,:,:]
    return out
