"""
Reads lesgo binary data and reformats into Python numpy format.
Author: Joel Bretheim
"""
import numpy as np
import re
from subprocess import check_output
from read_lesgo_bin import readmyfile
from os import getcwd, system

myDir = getcwd(); dirParts = myDir.split("/")
runName = dirParts[len(dirParts)-1]; print "This run's name: ", runName

dummy = check_output(["grep", 'nx, ny', "./lesgo_param.out"])
dummyStr = [int(s) for s in dummy.split() if s.isdigit()]
nx = dummyStr[0]; ny = dummyStr[1]; nz2 = dummyStr[2]; nz = dummyStr[3];
nz = nz - 1;
nz_ = nz2 - 1;
print nx, ny, nz, nz2, nz_

dummy = check_output(["grep", 'nproc', "./lesgo_param.out"])
dummyStr = [int(s) for s in dummy.split() if s.isdigit()]
nproc = dummyStr[0]
print nproc

#dummy = check_output(["grep", 'L_x', "./lesgo_param.out"])
#dummyStr = [float(s) for s in dummy.split() if s.isdigit()]
#a = re.findall(r"\d*([^]+)",dummyStr)
#print a

#kx_vec=[0,4,8]
#kx_vec=[0]
#kx_vec=[0,1,2]

avg       = 1;
snapshots = 0;  thisSnap = 5000;
fourier   = 0;

if avg:
    # neglect top point in each proc except last proc
    # remember, 0:8 does not grab the point at index 8, just indices 0-7
    vel  = np.zeros((3,nz,ny,nx))
    for i in range(0, nproc):
        fileName = './output/binary_vel_avg.dat.c' + str(i)
        filecontents = readmyfile(fileName)
        vel_i = np.reshape(filecontents, (3,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        vel[:,a:b,:,:] = vel_i[:,0:nz_,:,:]

    vel2  = np.zeros((6,nz,ny,nx))
    for i in range(0, nproc):
        fileName = './output/binary_vel2_avg.dat.c' + str(i)
        filecontents = readmyfile(fileName)
        vel_i = np.reshape(filecontents, (6,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        vel2[:,a:b,:,:] = vel_i[:,0:nz_,:,:]

    rs = np.zeros((6,nz,ny,nx))
    for i in range(0, nproc):
        fileName = './output/binary_rs.dat.c' + str(i)
        filecontents = readmyfile(fileName)
        rs_i = np.reshape(filecontents, (6,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        rs[:,a:b,:,:] = rs_i[:,0:nz_,:,:]

    tau = np.zeros((6,nz,ny,nx))
    for i in range(0, nproc):
        fileName = './output/binary_tau_avg.dat.c' + str(i)
        filecontents = readmyfile(fileName)
        tau_i = np.reshape(filecontents, (6,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        tau[:,a:b,:,:] = tau_i[:,0:nz_,:,:]

snap  = np.zeros((3,nz,ny,nx))
if snapshots:
    for i in range(0, nproc):
        fileName = './output/binary_vel.'+str(thisSnap)+'.dat.c'+str(i)
        filecontents = readmyfile(fileName)
        snap_i = np.reshape(filecontents, (3,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        snap[:,a:b,:,:] = snap_i[:,0:nz_,:,:]

if fourier:
    # complex-valued arrays corresponding to the real-valued arrays
    velc  = np.zeros((3,nz,ny,nx), dtype=complex)
    vel2c = np.zeros((6,nz,ny,nx), dtype=complex)
    rsc   = np.zeros((6,nz,ny,nx), dtype=complex)
    tauc  = np.zeros((6,nz,ny,nx), dtype=complex)
    snapc = np.zeros((3,nz,ny,nx), dtype=complex)

    # re-arrange real-valued arrays into complex-valued arrays
    for i in range(0,nx/2):
        b = 2*i+1;   a = b-1;   e = nx-i;  print 'i,a,b: ', i,a,b,e
        velc[:,:,:,i]  =  vel[:,:,:,a] + 1j *  vel[:,:,:,b]
        vel2c[:,:,:,i] = vel2[:,:,:,a] + 1j * vel2[:,:,:,b]
        rsc[:,:,:,i]   =   rs[:,:,:,a] + 1j *   rs[:,:,:,b]
        snapc[:,:,:,i] = snap[:,:,:,a] + 1j * snap[:,:,:,b]
        tauc[:,:,:,i]  =  tau[:,:,:,a] + 1j *  tau[:,:,:,b]
        if i > 0:
            velc[:,:,:,e] = np.conj(velc[:,:,:,i])
            vel2c[:,:,:,e] = np.conj(vel2c[:,:,:,i])
            rsc[:,:,:,e] = np.conj(rsc[:,:,:,i])

    for v in range(0,3):
        for i in range(0,nx/2):
            for k in range(0,nz):
                snapc[v,k,:,i] = np.fft.ifft(snapc[v,k,:,i]) * ny
                tauc[v,k,:,i] = np.fft.ifft(tauc[v,k,:,i]) * ny
                tauc[v+3,k,:,i] = np.fft.ifft(tauc[v+3,k,:,i]) * ny

    for i in range(1,nx/2):
        e = nx-i;  #print 'i,e : ', i,e
        snapc[:,:,:,e] = np.conj(snapc[:,:,:,i])
        tauc[:,:,:,e] = np.conj(tauc[:,:,:,i])

    velci = np.zeros((3,nz,ny,nx),dtype=complex)
    vel2ci = np.zeros((6,nz,ny,nx),dtype=complex)
    rsci = np.zeros((6,nz,ny,nx),dtype=complex)
    snapci = np.zeros((3,nz,ny,nx),dtype=complex)
    tauci = np.zeros((6,nz,ny,nx),dtype=complex)
    for v in range(0,3):
        for k in range(0,nz):
            for j in range(0,ny):
                velci[v,k,j,:] = np.fft.ifft(velc[v,k,j,:]) * nx
                vel2ci[v,k,j,:] = np.fft.ifft(vel2c[v,k,j,:]) * nx
                vel2ci[v+3,k,j,:] = np.fft.ifft(vel2c[v+3,k,j,:]) * nx
                rsci[v,k,j,:] = np.fft.ifft(rsc[v,k,j,:]) * nx
                rsci[v+3,k,j,:] = np.fft.ifft(rsc[v+3,k,j,:]) * nx
                snapci[v,k,j,:] = np.fft.ifft(snapc[v,k,j,:]) * nx
                tauci[v,k,j,:] = np.fft.ifft(tauc[v,k,j,:]) * nx
                tauci[v+3,k,j,:] = np.fft.ifft(tauc[v+3,k,j,:]) * nx

    # clean up and rename
    vel = np.real(velci)
    vel2 = np.real(vel2ci)
    rs = np.real(rsci)
    snap = np.real(snapci)
    tau = np.real(tauci)

else:
    print ">>>> Not in Fourier mode"

# rename data
u    =  vel[0,:,:,:]
v    =  vel[1,:,:,:]
w    =  vel[2,:,:,:]
uu   = vel2[0,:,:,:] 
vv   = vel2[1,:,:,:] 
ww   = vel2[2,:,:,:]
uw   = vel2[3,:,:,:] 
vw   = vel2[4,:,:,:] 
uv   = vel2[5,:,:,:]
rs11 =   rs[0,:,:,:] 
rs22 =   rs[1,:,:,:] 
rs33 =   rs[2,:,:,:]
rs13 =   rs[3,:,:,:] 
rs23 =   rs[4,:,:,:] 
rs12 =   rs[5,:,:,:]
txx  =  tau[0,:,:,:] 
txy  =  tau[1,:,:,:] 
tyy  =  tau[2,:,:,:]
txz  =  tau[3,:,:,:] 
tyz  =  tau[4,:,:,:] 
tzz  =  tau[5,:,:,:]

# compute horizontal averages
uXMean = np.mean(u, axis=2)      # x-average
uMean = np.mean(uXMean, axis=1)  # x- and y-averaged

txzMean = np.mean(txz, axis=2)      # x-averaging
txzMean = np.mean(txzMean, axis=1)  # y-averaging

rs11Mean = np.mean(rs11, axis=2)
rs11Mean = np.mean(rs11Mean, axis=1)
rs22Mean = np.mean(rs22, axis=2)
rs22Mean = np.mean(rs22Mean, axis=1)
rs33Mean = np.mean(rs33, axis=2)
rs33Mean = np.mean(rs33Mean, axis=1)
rs13Mean = np.mean(rs13, axis=2)
rs13Mean = np.mean(rs13Mean, axis=1)
rs23Mean = np.mean(rs23, axis=2)
rs23Mean = np.mean(rs23Mean, axis=1)
rs12Mean = np.mean(rs12, axis=2)
rs12Mean = np.mean(rs12Mean, axis=1)

uuMean = np.mean(uu, axis=2)
uuMean = np.mean(uuMean, axis=1)
vvMean = np.mean(vv, axis=2)
vvMean = np.mean(vvMean, axis=1)
wwMean = np.mean(ww, axis=2)
wwMean = np.mean(wwMean, axis=1)
uwMean = np.mean(uw, axis=2)
uwMean = np.mean(uwMean, axis=1)
vwMean = np.mean(vw, axis=2)
vwMean = np.mean(vwMean, axis=1)
uvMean = np.mean(uv, axis=2)
uvMean = np.mean(uvMean, axis=1)

datdir = 'pyData/'
system('mkdir ' + datdir)
np.save(datdir+'uXMean', uXMean)
np.save(datdir+'uMean', uMean)
np.save(datdir+'txzMean', txzMean)
np.save(datdir+'rs11Mean', rs11Mean)
np.save(datdir+'rs22Mean', rs22Mean)
np.save(datdir+'rs33Mean', rs33Mean)
np.save(datdir+'rs13Mean', rs13Mean)
np.save(datdir+'rs23Mean', rs23Mean)
np.save(datdir+'rs12Mean', rs12Mean)
np.save(datdir+'uuMean', uuMean)
np.save(datdir+'vvMean', vvMean)
np.save(datdir+'wwMean', wwMean)
np.save(datdir+'uwMean', uwMean)
np.save(datdir+'vwMean', vwMean)
np.save(datdir+'uvMean', uvMean)

np.save(datdir+'u', u)
np.save(datdir+'v', v)
np.save(datdir+'w', w)
np.save(datdir+'uu', uu)
np.save(datdir+'vv', vv)
np.save(datdir+'ww', ww)
np.save(datdir+'uw', uw)
np.save(datdir+'vw', vw)
np.save(datdir+'uv', uv)
np.save(datdir+'rs11', rs11)
np.save(datdir+'rs22', rs22)
np.save(datdir+'rs33', rs33)
np.save(datdir+'rs13', rs13)
np.save(datdir+'rs23', rs23)
np.save(datdir+'rs12', rs12)
np.save(datdir+'txx', txx)
np.save(datdir+'tyy', tyy)
np.save(datdir+'tzz', tzz)
np.save(datdir+'txz', txz)
np.save(datdir+'txy', txy)
np.save(datdir+'tyz', tyz)

np.save(datdir+'snap', snap)
