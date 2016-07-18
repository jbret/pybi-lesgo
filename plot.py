"""
Python plotting for lesgo binary data.
Author: Joel Bretheim
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('text', usetex = True)
import matplotlib.pyplot as plt
from matplotlib import ticker
import re
from subprocess import check_output
from read_lesgo_bin import readmyfile
from os import getcwd, system

RNL_branch = 1;    devel_branch = 0;

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
nz_ = nz2 - 1;
print "nx  =", nx
print "ny  =", ny
print "nz  =", nz
print "nz2 =", nz2
print "nz_ =", nz_

Lx = 2*np.pi;
Ly = 2*np.pi;
Lz = 1.0;

vel_avg_plot  = 1;
uXMean_plot   = 1;
tau_plot      = 1;
spanSpec_plot = 0;
sp1dky_plot   = 1;
sp1dkx_plot   = 1;
sp2d_plot     = 1;
rs_plot       = 1;
vel2_rs_plot  = 1;
snap_plot_xy  = 0;
snap_plot_yz  = 0;
snap_plot     = 0;  thisSnap = 5000;  # on uv-grid

mkm = 0;   # reference DNS data from MKM 1999

z = np.linspace(0, Lz, nz, endpoint=True)
y = np.linspace(0, Ly, ny, endpoint=False)
x = np.linspace(0, Lx, nx, endpoint=False)

datdir = 'data-npy/'
figdir1 = 'figs/'
figdir2 = figdir1 + 'otherFormats/'
figdir3 = figdir1 + 'highDPI/'
system('mkdir ' + figdir1)
system('mkdir ' + figdir2)
system('mkdir ' + figdir3)

def mySaveFig(figName, hiRes):
    print '>>> Saving ' + figName
    plt.savefig(figdir1 + figName + runName + '.png')
    plt.savefig(figdir2 + figName + runName + '.pdf')
    plt.savefig(figdir2 + figName + runName + '.eps')
    plt.savefig(figdir2 + figName + runName + '.jpg')
    if hiRes:
        myDPI = 600; dpiName = 'dpi'+str(myDPI)+'-';
        plt.savefig(figdir3+ dpiName + figName + runName + '.png', dpi=myDPI)
        plt.savefig(figdir3+ dpiName + figName + runName + '.pdf', dpi=myDPI)
        plt.savefig(figdir3+ dpiName + figName + runName + '.eps', dpi=myDPI)
        plt.savefig(figdir3+ dpiName + figName + runName + '.jpg', dpi=myDPI)

plt.close("all")
# begin plotting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if vel_avg_plot:
    uMean = np.load(datdir+'uMean.npy')
    fig = plt.figure()
    plt.semilogx(z, uMean, 'o')
    #plt.semilogx(z/(1./180), uMean, 'o')
    plt.semilogx(z, 1/0.4*np.log(z/.0001), '-k', label=r'$1/\kappa \ \mathrm{log}(z/z_{0})$')
    #plt.semilogx(z/(1./180), 1/0.41*np.log(z/(1./180))+5.0, '-k', label=r'$1/\kappa \ \mathrm{log}(z^{+})+B$')
    plt.xlabel('$ z / H $', fontsize=18);
    #plt.xlabel('$z^{+} $', fontsize=18);
    plt.ylabel('$[\ u / u_*]$', fontsize=18); 
    plt.xlim([.02 ,1.1])
    #plt.xlim([1, 1000])
    plt.text(.4,3,r'$ \kappa = 0.4,\ z_{0} = 10^{-4} $', fontsize=14)

    if mkm:
        yp180 = np.load('yp180.npy')
        Umean180 = np.load('Umean180.npy')
        plt.semilogx(yp180, Umean180,'o',label=r'$MKM99, DNS, Re_{\tau} = 180 $')

    plt.text(110,2.5,r'$ \kappa = 0.41,\ B = 5.0 $', fontsize=14)
    plt.legend(loc='lower right', fontsize=14)
    plt.tight_layout()

    mySaveFig('mvp_', 0)

if uXMean_plot:
    scale = 3.0;
    fig = plt.figure(figsize=(scale*Ly,scale*Lz))
    Y, Z = np.meshgrid(y, z)
    uXMean = np.load(datdir+'uXMean.npy')
    cs = plt.contourf(Y, Z, uXMean[:,:], vmin=0, vmax=17)
    cbar = plt.colorbar()
    plt.xlabel('$ y / H $', fontsize=18); plt.ylabel('$ z / H $', fontsize=18);
    #plt.suptitle('Streamwise velocity contours', fontsize = 16)
    # now make a circle with no fill, which is good for hilighting key results
    circle1=plt.Circle((0.261799, 0.1),.05,color='k',fill=False)
    circle2=plt.Circle((0.785398, 0.1),.05,color='k',fill=False)
    circle3=plt.Circle((1.309, 0.1),.05,color='k',fill=False)
    circle4=plt.Circle((1.8326, 0.1),.05,color='k',fill=False)
    circle5=plt.Circle((2.35619, 0.1),.05,color='k',fill=False)
    circle6=plt.Circle((2.87979, 0.1),.05,color='k',fill=False)
    ax = plt.gca()
    #ax.add_artist(circle1)
    #ax.add_artist(circle2)
    #ax.add_artist(circle3)
    #ax.add_artist(circle4)
    #ax.add_artist(circle5)
    #ax.add_artist(circle6)
    plt.tight_layout()
    mySaveFig('uXmean_', 0)

if tau_plot:
    rs13Mean = np.load(datdir+'rs13Mean.npy')
    txzMean = np.load(datdir+'txzMean.npy')
    fig = plt.figure()
    plt.plot(-1*rs13Mean, z, '-o', color='g', label = r'$[ -u^{\prime} w^{\prime}]$')
    plt.plot(-1*txzMean, z, '-o', color='r', label = r'$ [ -\tau_{xz} ] $')
    plt.plot(-1*(rs13Mean + txzMean), z, '-s', color='k', markeredgewidth=1, 
             markerfacecolor="None", label = r'$ \mathrm{sum} $')
    line = np.linspace(0,1,1000)
    plt.plot(line, 1+-1.0*line, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.15, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.05, '--', color='k')
    plt.xlabel(r'$ \mathrm{Stress} $', fontsize=18); plt.ylabel(r'$ z / H $', fontsize=18)
    plt.legend()
    plt.tight_layout()
    mySaveFig('tau_', 0)

if spanSpec_plot:
    sp1dky_uu = np.load(datdir+'sp1dky_uu.npy')
    sp1dky_uw = np.load(datdir+'sp1dky_uw.npy')
    sp1dky_vv = np.load(datdir+'sp1dky_vv.npy')
    sp1dky_ww = np.load(datdir+'sp1dky_ww.npy')

    spec11 = np.mean(sp1dky_uu[:,:,:], axis=2)
    spec13 = np.mean(sp1dky_uw[:,:,:], axis=2)
    spec22 = np.mean(sp1dky_vv[:,:,:], axis=2)
    spec33 = np.mean(sp1dky_ww[:,:,:], axis=2)

    #Euu180_0 = np.load('Euu180_0.npy')
    #Euu180_1 = np.load('Euu180_1.npy')
    #Euu180_2 = np.load('Euu180_2.npy')
    #Euu180_3 = np.load('Euu180_3.npy')
    #Evv180_0 = np.load('Evv180_0.npy')
    #Evv180_1 = np.load('Evv180_1.npy')
    #Evv180_2 = np.load('Evv180_2.npy')
    #Evv180_3 = np.load('Evv180_3.npy')
    #Eww180_0 = np.load('Eww180_0.npy')
    #Eww180_1 = np.load('Eww180_1.npy')
    #Eww180_2 = np.load('Eww180_2.npy')
    #Eww180_3 = np.load('Eww180_3.npy')

    heights = [ 5, 19, 30, 98 ];  numH = np.size(heights)
    ky = np.arange(0,ny/2)
    fig = plt.figure(figsize=(10,8))
    for j in range(0,numH):
        plt.subplot(1,numH,j+1)
        k = heights[j]
        plt.loglog(ky[1:], spec11[k,1:ny/2],'o', label=r'$E_{uu}$')
        plt.loglog(ky[1:], spec13[k,1:ny/2],'o', label=r'$E_{uw}$')
        plt.loglog(ky[1:], spec22[k,1:ny/2],'o', label=r'$E_{vv}$')
        plt.loglog(ky[1:], spec33[k,1:ny/2],'o', label=r'$E_{ww}$')
        if mkm:
            plt.loglog(yp180[1:], Euu180_0,'o', label=r'$DNS, E_{uu}$')
            plt.loglog(yp180[1:], Eww180_0,'o', label=r'$DNS, E_{vv}$')
            plt.loglog(yp180[1:], Evv180_0,'o', label=r'$DNS, E_{ww}$')
        
        plt.title(r'$ z^{+} = $'+str(heights[j]));
        plt.xlabel(r'$ k_y $'); plt.tight_layout()
        plt.legend(loc='lower left')
        #plt.ylim([])
    mySaveFig('spanSpec', 0)

if sp1dky_plot:
    sp1dky_uu = np.load(datdir+'sp1dky_uu.npy')
    sp1dky_uw = np.load(datdir+'sp1dky_uw.npy')
    sp1dky_vv = np.load(datdir+'sp1dky_vv.npy')
    sp1dky_ww = np.load(datdir+'sp1dky_ww.npy')

    e11 = np.mean(sp1dky_uu[:,:,:], axis=2)
    e13 = np.mean(sp1dky_uw[:,:,:], axis=2)
    e22 = np.mean(sp1dky_vv[:,:,:], axis=2)
    e33 = np.mean(sp1dky_ww[:,:,:], axis=2)

    ky = np.arange(0,ny/2)
    #lamY = Ly / ky
    lamY = ky

    #levels=[0.0,0.05,0.1,0.15,0.2,0.25,0.3]
    kyE11 = ky * e11[:,0:ny/2]
    kyE13 = ky * e13[:,0:ny/2]
    kyE22 = ky * e22[:,0:ny/2]
    kyE33 = ky * e33[:,0:ny/2]

    LAMY, Z = np.meshgrid(lamY[1:], z[1:])

    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(2,2,1)
    cs = plt.contourf(LAMY, Z, kyE11[1:,1:]);
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$k_y E_{uu} / u_{\tau}^2 $')
    plt.tight_layout()
     
    ax = fig.add_subplot(2,2,2)
    cs = plt.contourf(LAMY, Z, kyE11[1:,1:]);  
    #ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$k_y E_{uu} / u_{\tau}^2 $')
    plt.tight_layout()
    
    ax = fig.add_subplot(2,2,3)
    cs = plt.contourf(LAMY, Z, e11[1:,1:ny/2]);  
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$E_{uu} / u_{\tau}^2 $')
    plt.tight_layout()

    ax = fig.add_subplot(2,2,4)
    cs = plt.contourf(LAMY, Z, e11[1:,1:ny/2]);  
    #ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$E_{uu} / u_{\tau}^2 $')
    plt.tight_layout()

    mySaveFig('sp1dky_', 0)

    
if sp1dkx_plot:
    sp1dkx_uu = np.load(datdir+'sp1dkx_uu.npy')
    sp1dkx_uw = np.load(datdir+'sp1dkx_uw.npy')
    sp1dkx_vv = np.load(datdir+'sp1dkx_vv.npy')
    sp1dkx_ww = np.load(datdir+'sp1dkx_ww.npy')

    e11 = np.mean(sp1dkx_uu[:,:,:], axis=1)
    e13 = np.mean(sp1dkx_uw[:,:,:], axis=1)
    e22 = np.mean(sp1dkx_vv[:,:,:], axis=1)
    e33 = np.mean(sp1dkx_ww[:,:,:], axis=1)

    kx = np.arange(0,nx/2)
    #lamX = Lx / kx
    lamX = kx

    #levels=[0.0,0.05,0.1,0.15,0.2,0.25,0.3]
    kxE11 = kx * e11[:,0:nx/2]
    kxE13 = kx * e13[:,0:nx/2]
    kxE22 = kx * e22[:,0:nx/2]
    kxE33 = kx * e33[:,0:nx/2]

    LAMX, Z = np.meshgrid(lamX[1:], z[1:])

    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(2,2,1)
    cs = plt.contourf(LAMX, Z, kxE11[1:,1:]);
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_x / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$k_x E_{uu} / u_{\tau}^2 $')
    plt.tight_layout()
     
    ax = fig.add_subplot(2,2,2)
    cs = plt.contourf(LAMX, Z, kxE11[1:,1:]);  
    #ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_x / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$k_x E_{uu} / u_{\tau}^2 $')
    plt.tight_layout()
    
    ax = fig.add_subplot(2,2,3)
    cs = plt.contourf(LAMX, Z, e11[1:,1:nx/2]);  
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_x / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$E_{uu} / u_{\tau}^2 $')
    plt.tight_layout()

    ax = fig.add_subplot(2,2,4)
    cs = plt.contourf(LAMX, Z, e11[1:,1:nx/2]);  
    #ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_x / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$E_{uu} / u_{\tau}^2 $')
    plt.tight_layout()

    mySaveFig('sp1dkx_', 0)
    
if sp2d_plot:
    sp2d_uu = np.load(datdir+'sp2d_uu.npy')
    sp2d_uw = np.load(datdir+'sp2d_uw.npy')
    sp2d_vv = np.load(datdir+'sp2d_vv.npy')
    sp2d_ww = np.load(datdir+'sp2d_ww.npy')

    suu = sp2d_uu[:, 0:ny/2, 0:nx/2]
    suu_sum = np.sum(suu, axis=2)  # sum over all kx

    kx = np.arange(0,nx/2)
    ky = np.arange(0,ny/2)
    #lamY = Ly / ky
    lamY = ky

    #levels=[0.0,0.05,0.1,0.15,0.2,0.25,0.3]
    #sp11sum1 = np.sum(sp2d_uu[:,:,1:],axis=2)  # sum over kx
    #kyE11 = ky * sp11sum1[:,0:ny/2]

    #sp13sum1 = np.sum(sp2d_uw[:,:,1:],axis=2)  # sum over kx
    #kyE13 = ky * sp13sum1[:,0:ny/2]

    #sp22sum1 = np.sum(sp2d_vv[:,:,1:],axis=2)  # sum over kx
    #kyE22 = ky * sp22sum1[:,0:ny/2]

    #sp33sum1 = np.sum(sp2d_ww[:,:,1:],axis=2)  # sum over kx
    #kyE33 = ky * sp33sum1[:,0:ny/2]

    LAMY, Z = np.meshgrid(lamY[1:], z[1:])
    KX, KY = np.meshgrid(kx[1:], ky[1:])

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(2,2,1)
    #cs = plt.contourf(LAMY, Z, kyE11[1:,1:]);
    cs = plt.contourf(LAMY, Z, suu_sum[1:,1:]);
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$E_{uu} all k_x$')
    plt.tight_layout()

    ax = fig.add_subplot(2,2,2)
    #cs = plt.contourf(LAMY, Z, kyE13[1:,1:]); 
    cs = plt.contourf(LAMY, Z, suu[1:,1:,0]); 
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$E_{uu}, k_x=0$')
    plt.tight_layout()
     
    ax = fig.add_subplot(2,2,3)
    #cs = plt.contourf(LAMY, Z, kyE22[1:,1:]);  
    cs = plt.contourf(LAMY, Z, suu[1:,1:,9]);
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$E_{uu}, k_x=8$')
    plt.tight_layout()
    
    ax = fig.add_subplot(2,2,4)
    #cs = plt.contourf(LAMY, Z, kyE33[1:,1:]);
    cs = plt.contourf(LAMY, Z, suu[1:,1:,17]);  
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$E_{uu}, k_x=16$')
    plt.tight_layout()

    mySaveFig('sp2d_vert', 0)

    # horizontal >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(2,2,1)
    cs = plt.contourf(KY[1:20,1:5], KX[1:20,1:5], suu[4,1:20,1:5]);
    #ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_x / \delta $', fontsize=18);
    plt.ylabel(r'$ k_y / \delta $', fontsize=18); 
    plt.title(r'$E_{uu}$')
    plt.tight_layout()

    ax = fig.add_subplot(2,2,2)
    cs = plt.contourf(KY, KX, suu[8,1:,1:]);
    #ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_x / \delta $', fontsize=18);
    plt.ylabel(r'$ k_y / \delta $', fontsize=18); 
    plt.title(r'$E_{uu}$')
    plt.tight_layout()
     
    ax = fig.add_subplot(2,2,3)
    cs = plt.contourf(KY, KX, suu[20,1:,1:]);
    #ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_x / \delta $', fontsize=18);
    plt.ylabel(r'$ k_y / \delta $', fontsize=18); 
    plt.title(r'$E_{uu}$')
    plt.tight_layout()
    
    ax = fig.add_subplot(2,2,4)
    cs = plt.contourf(KY, KX, suu[44,1:,1:]);
    #ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ k_x / \delta $', fontsize=18);
    plt.ylabel(r'$ k_y / \delta $', fontsize=18); 
    plt.title(r'$E_{uu}$')
    plt.tight_layout()

    mySaveFig('sp2d_horiz', 0)
    
    
if rs_plot:
    rs11Mean = np.load(datdir+'rs11Mean.npy')
    rs22Mean = np.load(datdir+'rs22Mean.npy')
    rs33Mean = np.load(datdir+'rs33Mean.npy')
    rs13Mean = np.load(datdir+'rs13Mean.npy')
    rs23Mean = np.load(datdir+'rs23Mean.npy')
    rs12Mean = np.load(datdir+'rs12Mean.npy')
    fig = plt.figure()
    #plt.subplot(2,1,1)
    L1=plt.plot(z, rs11Mean, 'o', color='b', label = r'$[ u^{\prime} u^{\prime}]$')
    L2=plt.plot(z, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    L3=plt.plot(z, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    L4=plt.plot(z, rs13Mean, 'o', color='c', label = r'$[ u^{\prime} w^{\prime}]$')
    L5=plt.plot(z, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    L6=plt.plot(z, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    plt.legend(loc='upper right',ncol=2)
    plt.xlabel(r'$ z / H $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.tight_layout()
    #plt.subplot(2,1,2)
    #Re = 180; zp = z/(1./Re)
    #L1=plt.plot(zp, rs11Mean, 'o', color='b', label = r'$[ u^{\prime} u^{\prime}]$')
    #L2=plt.plot(zp, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    #L3=plt.plot(zp, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    #L4=plt.plot(zp, rs13Mean, 'o', color='c', label = r'$[ u^{\prime} w^{\prime}]$')
    #L5=plt.plot(zp, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    #L6=plt.plot(zp, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    #plt.legend(loc='upper right',ncol=2)
    #plt.xlabel(r'$ z^{+} $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    #plt.tight_layout()
    mySaveFig('rs_', 0)

    fig = plt.figure()
    #plt.subplot(2,1,1)
    L2=plt.plot(z, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    L3=plt.plot(z, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    L4=plt.plot(z, -1.0*rs13Mean, 'o', color='c', label = r'$[ -u^{\prime} w^{\prime}]$')
    L5=plt.plot(z, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    L6=plt.plot(z, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    plt.legend(loc='upper right',ncol=2)
    plt.xlabel(r'$ z / H $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.tight_layout()
    #plt.subplot(2,1,2)
    #Re = 180; zp = z/(1./Re)
    #L2=plt.plot(zp, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    #L3=plt.plot(zp, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    #L4=plt.plot(zp, -1.0*rs13Mean, 'o', color='c', label = r'$[ -u^{\prime} w^{\prime}]$')
    #L5=plt.plot(zp, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    #L6=plt.plot(zp, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    #plt.legend(loc='upper right',ncol=2)
    #plt.xlabel(r'$ z^{+} $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    #plt.tight_layout()
    mySaveFig('rs2_', 0)

    if vel2_rs_plot:
            fig = plt.figure()
            vvMean = np.load(datdir+'vvMean.npy')
            wwMean = np.load(datdir+'wwMean.npy')
            uwMean = np.load(datdir+'uwMean.npy')
            vwMean = np.load(datdir+'vwMean.npy')
            uvMean = np.load(datdir+'uvMean.npy')
            #L1=plt.plot(z, rs11Mean, 'o', color='b', label = r'$[ u^{\prime} u^{\prime}]$')
            L2=plt.plot(z, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
            L3=plt.plot(z, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
            L4=plt.plot(z, -1*rs13Mean, 'o', color='c', label = r'$[ u^{\prime} w^{\prime}]$')
            L5=plt.plot(z, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
            L6=plt.plot(z, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
            #L7=plt.plot(z, uuMean, color='b', label = r'$[ uu ]$')
            L8=plt.plot(z, vvMean, linewidth=3, color='g', label = r'$[ vv ]$')
            L9=plt.plot(z, wwMean, linewidth=3, color='r', label = r'$[ ww ]$')
            L10=plt.plot(z, -1*uwMean, linewidth=3, color='c', label = r'$[ uw ]$')
            L11=plt.plot(z, vwMean, linewidth=3, color='m', label = r'$[ vw ]$')
            L12=plt.plot(z, uvMean, linewidth=3, color='k', label = r'$[ uv ]$')
            plt.legend(loc='lower right',ncol=2)
            plt.xlabel(r'$ z / H $', fontsize=18); 
            plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
            plt.tight_layout()
            mySaveFig('vel2_rs_', 0)

            uuMean = np.load(datdir+'uuMean.npy')

            sp11_1d = np.load(datdir+'sp1dky_uu.npy')
            fig = plt.figure()
            plt.plot(z,rs11Mean,'^',label='rs11')
            plt.plot(z,uuMean,'s',label='uuMean')
            comp11 = np.mean(sp11_1d[:,:,:], axis=2)
            comp11_ = comp11[:,0:ny/2]
            comp11_[:,1:ny/2] = comp11_[:,1:ny/2] + comp11_[:,1:ny/2]
            comp11 = np.sum(comp11_[:,0:], axis=1)
            #comp11 = np.sum(comp11[:,0:], axis=1)
            plt.plot(z[1:], comp11[1:], 'o', label='spect')
            plt.legend()
            mySaveFig('comp11', 0)


            sp22_1d = np.load(datdir+'sp1dky_vv.npy')
            fig = plt.figure()
            plt.plot(z,rs22Mean,'^',label='rs22')
            plt.plot(z,vvMean,'s',label='vvMean')
            comp22 = np.mean(sp22_1d[:,:,:], axis=2)
            comp22_ = comp22[:,0:ny/2]
            comp22_[:,1:ny/2] = comp22_[:,1:ny/2] + comp22_[:,1:ny/2]
            comp22 = np.sum(comp22_[:,0:], axis=1)
            plt.plot(z[1:], comp22[1:], 'o', label='spect')
            plt.legend()
            mySaveFig('comp22', 0)
            
            sp33_1d = np.load(datdir+'sp1dky_ww.npy')
            fig = plt.figure()
            plt.plot(z,rs33Mean,'^',label='rs33')
            plt.plot(z,wwMean,'s',label='wwMean')
            comp33 = np.mean(sp33_1d[:,:,:], axis=2)
            comp33_ = comp33[:,0:ny/2]
            comp33_[:,1:ny/2] = comp33_[:,1:ny/2] + comp33_[:,1:ny/2]
            comp33 = np.sum(comp33_[:,0:], axis=1)
            plt.plot(z[1:], comp33[1:], 'o', label='spect')
            plt.legend()
            mySaveFig('comp33', 0)

if snap_plot_xy:
    snap = np.load(datdir+'snap.npy')
    X, Y = np.meshgrid(x, y)
    for k in range(2,3):
        scale = 3.0
        fig = plt.figure(figsize=(scale*Lx,scale*Ly))
        cs = plt.contourf(X, Y, snap[0,k,:,:])
        cbar = plt.colorbar()
        plt.xlabel(r'$ x / H $', fontsize=18); plt.ylabel(r'$ y / H $', fontsize=18); 
        plt.tight_layout()
        #fig.show()
        #plt.savefig(figdir+'xy_'+str(k)+'.png', dpi=100)

if snap_plot_yz:
    snap = np.load(datdir+'snap.npy')
    Y, Z = np.meshgrid(y, z)
    for i in range(2,3):
        scale = 3.0;
        fig = plt.figure(figsize=(scale*Ly,scale*Lz))
        cs = plt.contourf(Y, Z, snap[0,:,:,i]);  csName = 'yzCon_'
        #cs = plt.pcolor(Y, Z, snap[0,:,:,i]);  csName = 'yzCol_'
        cbar = plt.colorbar()
        plt.xlabel(r'$ y / H $', fontsize=18); plt.ylabel(r'$ z / H $', fontsize=18); 
        plt.tight_layout()
        mySaveFig(csName, 0)

        # print plt.gcf().canvas.get_supported_filetypes()
        # print plt.gcf().canvas.get_supported_filetypes_grouped()
        #
        # 'Scalable Vector Graphics': ['svg', 'svgz']
        # 'Portable Network Graphics': ['png']
        # 'Postscript': ['ps']
        # 'Joint Photographic Experts Group': ['jpeg', 'jpg']
        # 'Encapsulated Postscript': ['eps']
        # 'Tagged Image File Format': ['tif', 'tiff']
        # 'Raw RGBA bitmap': ['raw', 'rgba']
        # 'Portable Document Format': ['pdf']
        # 'Enhanced Metafile': ['emf']
