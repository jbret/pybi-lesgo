"""
Python plotting for lesgo binary data
Author: Joel Bretheim, jbretheim@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt
import re
from subprocess import check_output
from read_lesgo_bin import readmyfile
from os import getcwd

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

Lx = np.pi;
Ly = np.pi;
Lz = 1.0;
#kx_vec=[0,4,8]
kx_vec=[0]
#kx_vec=[0,1,2]

vel_avg_plot = 1;
rs_plot      = 1;
tau_plot     = 1;
snap_plot    = 0;  thisSnap = 5000;
snap_plot_yz = 0;
snap_plot_xy = 0;
fourier      = 1;

# neglect top point in each proc except last proc
# remember, 0:8 does not grab the point at index 8, just indices 0-7
vel  = np.zeros((3,nz,ny,nx))
# on w-grid
if vel_avg_plot:
    for i in range(0, nproc):
        fileName = './output/binary_vel_avg.dat.c' + str(i)
        filecontents = readmyfile(fileName)
        vel_i = np.reshape(filecontents, (3,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        vel[:,a:b,:,:] = vel_i[:,0:nz_,:,:]

rs = np.zeros((6,nz,ny,nx))
if rs_plot or tau_plot:
    for i in range(0, nproc):
        fileName = './output/binary_rs.dat.c' + str(i)
        filecontents = readmyfile(fileName)
        rs_i = np.reshape(filecontents, (6,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        rs[:,a:b,:,:] = rs_i[:,0:nz_,:,:]

tau = np.zeros((6,nz,ny,nx))
if tau_plot:
    for i in range(0, nproc):
        fileName = './output/binary_tau_avg.dat.c' + str(i)
        filecontents = readmyfile(fileName)
        tau_i = np.reshape(filecontents, (6,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        tau[:,a:b,:,:] = tau_i[:,0:nz_,:,:]

snap  = np.zeros((3,nz,ny,nx))
# on uv-grid
if snap_plot:
    for i in range(0, nproc):
        fileName = './output/binary_vel.'+str(thisSnap)+'.dat.c'+str(i)
        filecontents = readmyfile(fileName)
        snap_i = np.reshape(filecontents, (3,nz2,ny,nx))
        a = i*nz_;  b = (i+1)*nz_
        snap[:,a:b,:,:] = snap_i[:,0:nz_,:,:]

z = np.linspace(0, Lz, nz, endpoint=True)
y = np.linspace(0, Ly, ny, endpoint=False)
x = np.linspace(0, Lx, nx, endpoint=False)

if fourier:
    # complex-valued arrays corresponding to the real-valued arrays
    velc  = np.zeros((3,nz,ny,nx), dtype=complex)
    rsc   = np.zeros((6,nz,ny,nx), dtype=complex)
    tauc  = np.zeros((6,nz,ny,nx), dtype=complex)
    snapc = np.zeros((3,nz,ny,nx), dtype=complex)

    # re-arrange real-valued arrays into complex-valued arrays
    for i in range(0,nx/2):
        b = 2*i+1;   a = b-1;   e = nx-i;  print 'i,a,b: ', i,a,b,e
        velc[:,:,:,i]  =  vel[:,:,:,a] + 1j *  vel[:,:,:,b]
        rsc[:,:,:,i]   =   rs[:,:,:,a] + 1j *   rs[:,:,:,b]
        snapc[:,:,:,i] = snap[:,:,:,a] + 1j * snap[:,:,:,b]
        tauc[:,:,:,i]  =  tau[:,:,:,a] + 1j *  tau[:,:,:,b]
        if i > 0:
            velc[:,:,:,e] = np.conj(velc[:,:,:,i])
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
    rsci = np.zeros((6,nz,ny,nx),dtype=complex)
    snapci = np.zeros((3,nz,ny,nx),dtype=complex)
    tauci = np.zeros((6,nz,ny,nx),dtype=complex)
    for v in range(0,3):
        for k in range(0,nz):
            for j in range(0,ny):
                velci[v,k,j,:] = np.fft.ifft(velc[v,k,j,:]) * nx
                rsci[v,k,j,:] = np.fft.ifft(rsc[v,k,j,:]) * nx
                rsci[v+3,k,j,:] = np.fft.ifft(rsc[v+3,k,j,:]) * nx
                snapci[v,k,j,:] = np.fft.ifft(snapc[v,k,j,:]) * nx
                tauci[v,k,j,:] = np.fft.ifft(tauc[v,k,j,:]) * nx
                tauci[v+3,k,j,:] = np.fft.ifft(tauc[v+3,k,j,:]) * nx

    # clean up and rename
    vel = np.real(velci)
    rs = np.real(rsci)
    snap = np.real(snapci)
    tau = np.real(tauci)

else:
    print ">>>> Not in Fourier mode"

# rename data
u    = vel[0,:,:,:]
v    = vel[1,:,:,:]
w    = vel[2,:,:,:]
rs11 =  rs[0,:,:,:] 
rs22 =  rs[1,:,:,:] 
rs33 =  rs[2,:,:,:]
rs13 =  rs[3,:,:,:] 
rs23 =  rs[4,:,:,:] 
rs12 =  rs[5,:,:,:]
txx  = tau[0,:,:,:] 
txy  = tau[1,:,:,:] 
tyy  = tau[2,:,:,:]
txz  = tau[3,:,:,:] 
tyz  = tau[4,:,:,:] 
tzz  = tau[5,:,:,:]

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

plt.close("all")
# begin plotting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if vel_avg_plot:
    fig = plt.figure()
    plt.semilogx(z, uMean, 'o')
    plt.semilogx(z, 1/0.4*np.log(z/.0001), '-k', label=r'$1/\kappa \ \mathrm{log}(z/z_{0})$')
    plt.xlabel('$z / H$', fontsize=18); plt.ylabel('$[\ u / u_*]$', fontsize=18); 
    plt.xlim([.02 ,1.1])
    plt.text(.4,3,r'$ \kappa = 0.4,\ z_{0} = 10^{-4} $', fontsize=14)
    plt.legend(loc='lower right', fontsize=14)
    plt.tight_layout()
    plt.savefig('mvp_' + runName + '.png')
    fig.show()

    scale = 3.0;
    fig = plt.figure(figsize=(scale*Ly,scale*Lz))
    Y, Z = np.meshgrid(y, z)
    cs = plt.contourf(Y, Z, uXMean[:,:])
    cbar = plt.colorbar()
    plt.xlabel('$ y / H $', fontsize=18); plt.ylabel('$ z / H $', fontsize=18);
    #plt.suptitle('Streamwise velocity contours', fontsize = 16)
    # now make a circle with no fill, which is good for hilighting key results
    circle1=plt.Circle((y[1*ny/5],0.1),.05,color='k',fill=False)
    circle2=plt.Circle((y[2*ny/5],0.1),.05,color='k',fill=False)
    circle3=plt.Circle((y[3*ny/5],0.1),.05,color='k',fill=False)
    circle4=plt.Circle((y[4*ny/5],0.1),.05,color='k',fill=False)
    ax = plt.gca()
    ax.add_artist(circle1)
    ax.add_artist(circle2)
    ax.add_artist(circle3)
    ax.add_artist(circle4)
    plt.tight_layout()
    plt.savefig('uXmean_' + runName + '.png')  #,dpi = 300
    fig.show()

if tau_plot:
    fig = plt.figure()
    plt.plot(-1*rs13Mean, z, '-o', color='g', label = r'$[ u^{\prime} w^{\prime}]$')
    plt.plot(-1*txzMean, z, '-o', color='r', label = r'$ [ \tau_{xz} ] $')
    plt.plot(-1*(rs13Mean + txzMean), z, '-s', color='k', markeredgewidth=1, markerfacecolor="None", label = r'$ \mathrm{sum} $')
    line = np.linspace(0,1,1000)
    plt.plot(line, 1+-1.0*line, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.15, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.05, '--', color='k')
    plt.xlabel(r'$ \mathrm{Stress} $', fontsize=18); plt.ylabel(r'$ z / H $', fontsize=18)
    plt.legend()
    plt.tight_layout()
    plt.savefig('tau_' + runName + '.png')
    fig.show()
    
if rs_plot:
    fig = plt.figure()
    L1=plt.plot(z, rs11Mean, 'o', color='b', label = r'$[ u^{\prime} u^{\prime}]$')
    L2=plt.plot(z, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    L3=plt.plot(z, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    L4=plt.plot(z, -1*rs13Mean, 'o', color='c', label = r'$[ u^{\prime} w^{\prime}]$')
    L5=plt.plot(z, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    L6=plt.plot(z, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    plt.legend(loc='upper right')
    plt.xlabel(r'$ z / H $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.tight_layout()
    plt.savefig('rs_' + runName + '.png')
    fig.show()

if snap_plot_xy:
    X, Y = np.meshgrid(x, y)
    for k in range(2,3):
        scale = 3.0
        fig = plt.figure(figsize=(scale*Lx,scale*Ly))
        cs = plt.contourf(X, Y, snap[0,k,:,:])
        cbar = plt.colorbar()
        plt.xlabel(r'$ x / H $', fontsize=18); plt.ylabel(r'$ y / H $', fontsize=18); 
        plt.tight_layout()
        fig.show()
        #plt.savefig('xy_'+str(k)+'.png', dpi=100)

if snap_plot_yz:
    Y, Z = np.meshgrid(y, z)
    for i in range(2,3):
        scale = 3.0;
        fig = plt.figure(figsize=(scale*Ly,scale*Lz))
        cs = plt.contourf(Y, Z, snap[0,:,:,i]);  csName = 'yzCon'
        #cs = plt.pcolor(Y, Z, snap[0,:,:,i]);  csName = 'yzCol'
        cbar = plt.colorbar()
        plt.xlabel(r'$ y / H $', fontsize=18); plt.ylabel(r'$ z / H $', fontsize=18); 
        plt.tight_layout()
        plt.savefig(csName + runName + '.png') #,dpi=300)
        fig.show()
        
        #plt.savefig('yzCont_'+str(i)+'.png', dpi=300)
        #plt.savefig('yzPcol_'+str(i)+'.png')
        #plt.savefig('yzPcol_'+str(i)+'.pdf', dpi=100)
        #plt.savefig('yzPcol_'+str(i)+'.eps', dpi=100)
        #print plt.gcf().canvas.get_supported_filetypes()
        #print plt.gcf().canvas.get_supported_filetypes_grouped()
        # 'Scalable Vector Graphics': ['svg', 'svgz']
        # 'Portable Network Graphics': ['png']
        # 'Postscript': ['ps']
        # 'Joint Photographic Experts Group': ['jpeg', 'jpg']
        # 'Encapsulated Postscript': ['eps']
        # 'Tagged Image File Format': ['tif', 'tiff']
        # 'Raw RGBA bitmap': ['raw', 'rgba']
        # 'Portable Document Format': ['pdf']
        # 'Enhanced Metafile': ['emf']

        

