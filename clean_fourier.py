"""
Reads lesgo's binary data output into a numpy array
Author: Joel Bretheim
"""
import numpy as np

def clean_fourier(u):
    myShape = np.shape(u)
    print 'shape: ',  myShape
    uc = np.zeros(myShape, dtype=complex) # (kx,y,z)
    
    # re-arrange real-valued arrays into complex-valued arrays
    for i in range(0,nx/2):
        b = 2*i+1;   a = b-1;   e = nx-i;  print 'i,a,b: ', i,a,b,e
        uc[:,:,:,i]  =  u[:,:,:,a] + 1j*u[:,:,:,b]
        if i > 0:
            uc[:,:,:,e] = np.conj( uc[:,:,:,i] )
  
    # go from kx space to x space
    uci = np.zeros(myShape, dtype=complex)
    for v in range(0,3):
        for k in range(0,nz):
            for j in range(0,ny):
                uci[v,k,j,:] = np.fft.ifft(uc[v,k,j,:]) * nx

    u_clean = np.real(uci)
    return u_clean
