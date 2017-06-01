"""
Reads lesgo binary data and reformats into Python numpy format.
Author: Joel Bretheim
"""
import numpy as np
import re
from subprocess import check_output
#from read_lesgo_bin import readmyfile
from os import getcwd, system
from assemble_field import assemble_field
from unfoldings import unfold

RNL_branch = 1;    devel_branch = 0;
avg        = 1;
snapshots  = 0;    thisSnap = 250300;
fourier    = 0;
spectra_jb = 0;

myDir = getcwd(); dirParts = myDir.split("/")
runName = dirParts[len(dirParts)-1]; print "This run's name: ", runName

if RNL_branch:
    lesgo_param_loc = "./lesgo_param.out"
elif devel_branch:
    lesgo_param_loc = "./output/lesgo_param.out"
else:
    print "Must specify location of lesgo_param.out"

dummy = check_output(["grep", 'nx, ny', lesgo_param_loc])
dummyStr = [int(s) for s in dummy.split() if s.isdigit()]
nx = dummyStr[0]; ny = dummyStr[1]; nz2 = dummyStr[2]; nz = dummyStr[3];
nz = nz - 1;
#nz_ = nz2 - 1;
print "nx  =", nx
print "ny  =", ny
print "nz  =", nz
print "nz2 =", nz2
#print "nz_ =", nz_

dummy = check_output(["grep", 'nproc', lesgo_param_loc])
dummyStr = [int(s) for s in dummy.split() if s.isdigit()]
nproc = dummyStr[0]
print "nproc =", nproc

print "nx*ny*nz2*8*1 =", nx*ny*nz2*8*1
print "nx*ny*nz2*8*3 =", nx*ny*nz2*8*3
print "nx*ny*nz2*8*6 =", nx*ny*nz2*8*6

if avg:
    vel    = assemble_field('./output/binary_vel_avg.dat.c',nproc,3,nz2,nz,ny,nx)
    vel2   = assemble_field('./output/binary_vel2_avg.dat.c',nproc,6,nz2,nz,ny,nx)
    rs     = assemble_field('./output/binary_rs.dat.c',nproc,6,nz2,nz,ny,nx)
    tau    = assemble_field('./output/binary_tau_avg.dat.c',nproc,6,nz2,nz,ny,nx)
    dp     = assemble_field('./output/binary_dp.dat.c',nproc,3,nz2,nz,ny,nx)
    f      = assemble_field('./output/binary_force_avg.dat.c',nproc,3,nz2,nz,ny,nx)
    sp2d   = assemble_field('./output/binary_sp2d.dat.c',nproc,6,nz2,nz,ny,nx)
    sp1dky = assemble_field('./output/binary_sp1dky.dat.c',nproc,6,nz2,nz,ny,nx)
    sp1dkx = assemble_field('./output/binary_sp1dkx.dat.c',nproc,6,nz2,nz,ny,nx)
    spvort = np.zeros((8,nz,ny,nx))
    nu_t = np.zeros((1,nz,ny,nx))
    if fourier == 0:
        spvort = assemble_field('./output/binary_spvort.dat.c',nproc,8,nz2,nz,ny,nx)
        nu_t   = assemble_field('./output/binary_nu_t.dat.c',nproc,1,nz2,nz,ny,nx)
    snap = np.zeros((3,nz,ny,nx))
    if snapshots:
        snap   = assemble_field('./output/binary_vel.'+str(thisSnap)+'.dat.c',
                 nproc,3,nz2,nz,ny,nx)
if fourier:
    vel  = unfold(vel,  1)   # kx --> x
    vel2 = unfold(vel2, 1)
    rs   = unfold(rs,   1)
    tau  = unfold(tau,  2)  # kx, ky --> x
    snap = unfold(snap, 2)  # not tested yet

if spectra_jb:
    # complex-valued arrays corresponding to the real-valued arrays
    # stored in (kx,ky) space
    sp2dc   = np.zeros((6,nz,ny,nx), dtype=complex) # (kx,ky,z)
    sp2dci  = np.zeros((6,nz,ny,nx), dtype=complex) # (kx,ky,z)
    # stored in (kx) space
    sp1dkxc   = np.zeros((6,nz,ny,nx), dtype=complex) # (kx, y,z)
    # only the last 4 entries are in kx space
    spvortc   = np.zeros((4,nz,ny,nx), dtype=complex) # (kx, y,z)

    # re-arrange real-valued arrays into complex-valued arrays
    for i in range(0,nx/2):
        b = 2*i+1;   a = b-1;   e = nx-i;  print 'i,a,b: ', i,a,b,e
        sp1dkxc[:,:,:,i]  =  sp1dkx[:,:,:,a] + 1j * sp1dkx[:,:,:,b]
        # only the fifth entry (index 4) is in kx space
        spvortc[0:4,:,:,i]  =  spvort[4:8,:,:,a] + 1j * spvort[4:8,:,:,b]
        sp2dc[:,:,:,i]    =  sp2d[:,:,:,a]   + 1j * sp2d[:,:,:,b]

    # go from kx,ky space to kx space
    #for v in range(0,6):
    #    for i in range(0,nx/2):
    #        for k in range(0,nz):
    #            sp2dci[v,k,:,i] = np.fft.ifft(sp2dc[v,k,:,i]) * ny
    
    #for i in range(1,nx/2):
    #    e = nx-i;  #print 'i,e : ', i,e
    #    sp2dci[:,:,:,e] = np.conj(sp2dci[:,:,:,i])

    # clean up and rename
    #snap = np.real(snapci)
    sp1dkx = np.real(sp1dkxc)
    spvort[4:8,:,:,:] = np.real(spvortc[0:4,:,:,:])
    sp2d   = np.real(sp2dc)
    
else:
    print ">>>> No spectra_jb"



# rename data
u    =  vel[0,:,:,:]
v    =  vel[1,:,:,:]
w    =  vel[2,:,:,:]
dpdx = dp[0,:,:,:]
dpdy = dp[1,:,:,:]
dpdz = dp[2,:,:,:]
fx = f[0,:,:,:]
fy = f[1,:,:,:]
fz = f[2,:,:,:]
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
sp2d_uu =   sp2d[0,:,:,:] 
sp2d_vv =   sp2d[1,:,:,:] 
sp2d_ww =   sp2d[2,:,:,:]
sp2d_uw =   sp2d[3,:,:,:] 
sp2d_vw =   sp2d[4,:,:,:] 
sp2d_uv =   sp2d[5,:,:,:]
sp1dky_uu =   sp1dky[0,:,:,:] 
sp1dky_vv =   sp1dky[1,:,:,:] 
sp1dky_ww =   sp1dky[2,:,:,:]
sp1dky_uw =   sp1dky[3,:,:,:] 
sp1dky_vw =   sp1dky[4,:,:,:] 
sp1dky_uv =   sp1dky[5,:,:,:]
sp1dkx_uu =   sp1dkx[0,:,:,:] 
sp1dkx_vv =   sp1dkx[1,:,:,:] 
sp1dkx_ww =   sp1dkx[2,:,:,:]
sp1dkx_uw =   sp1dkx[3,:,:,:] 
sp1dkx_vw =   sp1dkx[4,:,:,:] 
sp1dkx_uv =   sp1dkx[5,:,:,:]
spvort_vortx = spvort[0,:,:,:] 
spvort_vorty = spvort[1,:,:,:] 
spvort_vortz = spvort[2,:,:,:] 
spvort_vortp = spvort[3,:,:,:] 
spvort_vortsx = spvort[4,:,:,:] 
spvort_vortsy = spvort[5,:,:,:] 
spvort_vortsz = spvort[6,:,:,:] 
spvort_vorts = spvort[7,:,:,:] 
txx  =  tau[0,:,:,:] 
txy  =  tau[1,:,:,:] 
tyy  =  tau[2,:,:,:]
txz  =  tau[3,:,:,:] 
tyz  =  tau[4,:,:,:] 
tzz  =  tau[5,:,:,:]
nu_t  =  nu_t[0,:,:,:] 

# compute horizontal averages
uXMean = np.mean(u, axis=2)      # x-average
vXMean = np.mean(v, axis=2)      # x-average
wXMean = np.mean(w, axis=2)      # x-average
uMean = np.mean(uXMean, axis=1)  # x- and y-averaged
vMean = np.mean(vXMean, axis=1)  # x- and y-averaged
wMean = np.mean(wXMean, axis=1)  # x- and y-averaged

nu_tMean = np.mean(nu_t, axis=2) # x-average
nu_tMean = np.mean(nu_tMean, axis=1) # x- and y-averaged

txzMean = np.mean(txz, axis=2)      # x-averaging
txzMean = np.mean(txzMean, axis=1)  # y-averaging

dpdxMean = np.mean(dpdx, axis=2)      # x-averaging
dpdxMean = np.mean(dpdxMean, axis=1)  # y-averaging
dpdyMean = np.mean(dpdy, axis=2)      # x-averaging
dpdyMean = np.mean(dpdyMean, axis=1)  # y-averaging
dpdzMean = np.mean(dpdz, axis=2)      # x-averaging
dpdzMean = np.mean(dpdzMean, axis=1)  # y-averaging

fxMean = np.mean(fx, axis=2)      # x-averaging
fxMean = np.mean(fxMean, axis=1)  # y-averaging
fyMean = np.mean(fy, axis=2)      # x-averaging
fyMean = np.mean(fyMean, axis=1)  # y-averaging
fzMean = np.mean(fz, axis=2)      # x-averaging
fzMean = np.mean(fzMean, axis=1)  # y-averaging

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

datdir = 'data-npy/'
system('mkdir ' + datdir)
np.save(datdir+'uXMean', uXMean)
np.save(datdir+'uMean', uMean)
np.save(datdir+'vMean', vMean)
np.save(datdir+'wMean', wMean)
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

np.save(datdir+'fx', fx)
np.save(datdir+'fy', fy)
np.save(datdir+'fz', fz)

np.save(datdir+'fxMean', fxMean)
np.save(datdir+'fyMean', fyMean)
np.save(datdir+'fzMean', fzMean)
np.save(datdir+'dpdxMean', dpdxMean)
np.save(datdir+'dpdyMean', dpdyMean)
np.save(datdir+'dpdzMean', dpdzMean)

np.save(datdir+'dpdx', dpdx)
np.save(datdir+'dpdy', dpdy)
np.save(datdir+'dpdz', dpdz)

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
np.save(datdir+'sp2d_uu', sp2d_uu)
np.save(datdir+'sp2d_vv', sp2d_vv)
np.save(datdir+'sp2d_ww', sp2d_ww)
np.save(datdir+'sp2d_uw', sp2d_uw)
np.save(datdir+'sp2d_vw', sp2d_vw)
np.save(datdir+'sp2d_uv', sp2d_uv)
np.save(datdir+'sp1dky_uu', sp1dky_uu)
np.save(datdir+'sp1dky_vv', sp1dky_vv)
np.save(datdir+'sp1dky_ww', sp1dky_ww)
np.save(datdir+'sp1dky_uw', sp1dky_uw)
np.save(datdir+'sp1dky_vw', sp1dky_vw)
np.save(datdir+'sp1dky_uv', sp1dky_uv)
np.save(datdir+'sp1dkx_uu', sp1dkx_uu)
np.save(datdir+'sp1dkx_vv', sp1dkx_vv)
np.save(datdir+'sp1dkx_ww', sp1dkx_ww)
np.save(datdir+'sp1dkx_uw', sp1dkx_uw)
np.save(datdir+'sp1dkx_vw', sp1dkx_vw)
np.save(datdir+'sp1dkx_uv', sp1dkx_uv)
np.save(datdir+'spvort_vortx', spvort_vortx)
np.save(datdir+'spvort_vorty', spvort_vorty)
np.save(datdir+'spvort_vortz', spvort_vortz)
np.save(datdir+'spvort_vortp', spvort_vortp)
np.save(datdir+'spvort_vortsx', spvort_vortsx)
np.save(datdir+'spvort_vortsy', spvort_vortsy)
np.save(datdir+'spvort_vortsz', spvort_vortsz)
np.save(datdir+'spvort_vorts', spvort_vorts)
np.save(datdir+'txx', txx)
np.save(datdir+'tyy', tyy)
np.save(datdir+'tzz', tzz)
np.save(datdir+'txz', txz)
np.save(datdir+'txy', txy)
np.save(datdir+'tyz', tyz)

np.save(datdir+'nu_t', nu_t)
np.save(datdir+'nu_tMean', nu_tMean)

np.save(datdir+'snap', snap)
