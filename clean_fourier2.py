"""
Reads lesgo's binary data output into a numpy array
Author: Joel Bretheim
"""
import numpy as np

def clean_fourier2(u):
    myShape = np.shape(u); print 'shape: ', myShape
    num = myShape[0]; print 'num: ', num
    nz  = myShape[1]; print 'nz: ', nz
    ny  = myShape[2]; print 'ny: ', ny
    nx  = myShape[3]; print 'nx: ', nx

    uc = np.zeros(myShape, dtype=complex) # (kx,y,z)
    
    # re-arrange real-valued arrays into complex-valued arrays
    for i in range(0,nx/2):
        b = 2*i+1;   a = b-1;   e = nx-i;  print 'i,a,b: ', i,a,b,e
        uc[:,:,:,i]  =  u[:,:,:,a] + 1j*u[:,:,:,b]

    # go from kx,ky space to kx space
    for v in range(0,3): # should indeed be 3 and not num
        for i in range(0,nx/2):
            for k in range(0,nz):
                uc[v,k,:,i] = np.fft.ifft(uc[v,k,:,i]) * ny
                if num == 6:
                    uc[v+3,k,:,i] = np.fft.ifft(uc[v+3,k,:,i]) * ny

    for i in range(1,nx/2):
        e = nx-i;  #print 'i,e : ', i,e
        uc[:,:,:,e] = np.conj(uc[:,:,:,i])

    uci = np.zeros((num,nz,ny,nx),dtype=complex)

    # go from kx space to x space
    for v in range(0,3):
        for k in range(0,nz):
            for j in range(0,ny):
                uci[v,k,j,:] = np.fft.ifft(uc[v,k,j,:]) * nx
                if num == 6:
                    uci[v+3,k,j,:] = np.fft.ifft(uc[v+3,k,j,:]) * nx

    u_clean = np.real(uci)
    return u_clean

