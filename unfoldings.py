"""
Reads lesgo's binary data output into a numpy array
Author: Joel Bretheim
"""
import numpy as np

def getParams(u):
    myShape = np.shape(u);
    num = myShape[0]; nz  = myShape[1]; ny  = myShape[2]; nx  = myShape[3]; 
    print 'params (shape, num, nz, ny, nx): ', myShape, num, nz, ny, nx
    return myShape, num, nz, ny, nx

def kx2x(uc):
    myShape, num, nz, ny, nx = getParams(uc)
    # go from kx space to x space
    uci = np.zeros(myShape, dtype=complex)
    for v in range(0,3):  # should indeed be 3 and not num
        for k in range(0,nz):
            for j in range(0,ny):
                uci[v,k,j,:] = np.fft.ifft(uc[v,k,j,:]) * nx
                if num == 6:
                    uci[v+3,k,j,:] = np.fft.ifft(uc[v+3,k,j,:]) * nx

    return np.real(uci)

def kxky2kx(uc):
    myShape, num, nz, ny, nx = getParams(uc)
    # go from kx,ky space to kx space
    for v in range(0,3): # should indeed be 3 and not num
        for i in range(0,nx/2):
            for k in range(0,nz):
                uc[v,k,:,i] = np.fft.ifft(uc[v,k,:,i]) * ny
                if num == 6:
                    uc[v+3,k,:,i] = np.fft.ifft(uc[v+3,k,:,i]) * ny

    for i in range(1,nx/2):
        e = nx-i;
        uc[:,:,:,e] = np.conj(uc[:,:,:,i])
    return uc

def unfold(u, myType):
    myShape, num, nz, ny, nx = getParams(u)
    uc = np.zeros(myShape, dtype=complex) # (kx,y,z)
    
    # re-arrange real-valued arrays into complex-valued arrays
    for i in range(0,nx/2):
        b = 2*i+1;   a = b-1;   e = nx-i;
        uc[:,:,:,i]  =  u[:,:,:,a] + 1j*u[:,:,:,b]
        if myType ==1 and i > 0:
            uc[:,:,:,e] = np.conj( uc[:,:,:,i] )

    if myType == 2:
        # go from kx,ky space to kx space
        uc = kxky2kx(uc)
    
    out = kx2x(uc)
    return out

