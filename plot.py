"""
Python plotting for lesgo binary data.
Author: Joel Bretheim
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
#matplotlib.rc('text', usetex = True)   # may need to disable on some machines
import matplotlib.pyplot as plt
from matplotlib import ticker
import re
from subprocess import check_output
from read_lesgo_bin import readmyfile
from os import getcwd, system

RNL_branch = 1;    devel_branch = 0;

vel_avg_plot    = 0;
uXMean_plot     = 0;
tau_plot        = 0;
spanSpec_plot   = 0;
sp1dky_plot     = 1;  plot_wavelen = 0;  # by wavelength or wavenumber
sp1dkx_plot     = 1;
sp2d_plot_vert  = 0;  # WARNING: must test plot labels in LES case for sp2d plots (both vert and horiz)
sp2d_plot_horiz = 0;  localMax = 1  # if 0 then uses global max
rs_plot         = 0;
vel2_rs_plot    = 0;
snap_plot_xy    = 0;
snap_plot_yz    = 0;
snap_plot       = 0;  thisSnap = 5000;  # on uv-grid

mkm = 0;   # reference DNS data from MKM 1999

vert = 'z';  vertvel = 'w';
span = 'y';  spanvel = 'v';

compare_cases = 0;
c1dir = '../rnl-64-kx1/'

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

z = np.linspace(0, Lz, nz, endpoint=False)  # w-grid (no point at z = 1.0)
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
    plt.semilogx(z, uMean, 'or') #, label=r'$ k_x = 2 $')
    #plt.semilogx(z/(1./180), uMean, 'o')
    if compare_cases:
        c1_uMean = np.load(c1dir+datdir+'uMean.npy')
        c1_z = z
        plt.semilogx(c1_z, c1_uMean, 'or', label=r'$ k_x = 3 $')
    plt.semilogx(z, 1/0.4*np.log(z/.0001), '-k', 
                 label=r'$1/\kappa \ \mathrm{log}('+vert+'/'+vert+'_{0})$')
    #plt.semilogx(z/(1./180), 1/0.41*np.log(z/(1./180))+5.0, '-k', label=r'$1/\kappa \ \mathrm{log}(z^{+})+B$')
    plt.xlabel('$'+vert+' / H $', fontsize=18);
    #plt.xlabel('$'+vert+'^{+} $', fontsize=18);
    plt.ylabel('$[\ u / u_*]$', fontsize=18); 
    #plt.xlim([.02 ,1.1])
    #plt.xlim([1, 1000])
    plt.text(.3,3,r'$ \kappa = 0.4,\ '+vert+'_{0} = 10^{-4} $', fontsize=14)
    #plt.text(.06,1,r'$ \kappa = 0.4,\ '+vert+'_{0} = 10^{-4} $', fontsize=14)
    plt.xticks([10**(-2), 10**(-1), 10**0], fontsize = 16)
    plt.yticks([0,5,10,15,20,25], fontsize = 16)

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
    plt.xlabel('$'+span+' / H $', fontsize=18); 
    plt.ylabel('$'+vert+' / H $', fontsize=18);
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
    plt.plot(-1*rs13Mean, z, '-o', color='g', 
              label = '$[ -u^{\prime}'+vertvel+'^{\prime}]$')
    # note below the 'r' needed before label (in order to render the \tau)
    plt.plot(-1*txzMean, z, '-o', color='r', label = r'$ [ -\tau_{x'+vert+'} ] $')
    plt.plot(-1*(rs13Mean + txzMean), z, '-s', color='k', markeredgewidth=1, 
             markerfacecolor="None", label = r'$ \mathrm{sum} $')
    line = np.linspace(0,1,1000)
    plt.plot(line, 1+-1.0*line, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.15, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.05, '--', color='k')
    plt.xlabel('$ \mathrm{Stress} $', fontsize=18); 
    plt.ylabel('$'+vert+' / H $', fontsize=18)
    plt.legend()
    plt.tight_layout()
    mySaveFig('tau_', 0)

if spanSpec_plot:
    sp1dky_uu = np.load(datdir+'sp1dky_uu.npy')
    #sp1dky_uw = np.load(datdir+'sp1dky_uw.npy')
    sp1dky_vv = np.load(datdir+'sp1dky_vv.npy')
    sp1dky_ww = np.load(datdir+'sp1dky_ww.npy')

    spec11 = np.mean(sp1dky_uu[:,:,:], axis=2)
    #spec13 = np.mean(sp1dky_uw[:,:,:], axis=2)
    spec22 = np.mean(sp1dky_vv[:,:,:], axis=2)
    spec33 = np.mean(sp1dky_ww[:,:,:], axis=2)

    if mkm:
        Euu180_0z = np.load('Euu180_0z.npy')
        Euu180_1z = np.load('Euu180_1z.npy')
        Euu180_2z = np.load('Euu180_2z.npy')
        Evv180_0z = np.load('Evv180_0z.npy')
        Evv180_1z = np.load('Evv180_1z.npy')
        Evv180_2z = np.load('Evv180_2z.npy')
        Eww180_0z = np.load('Eww180_0z.npy')
        Eww180_1z = np.load('Eww180_1z.npy')
        Eww180_2z = np.load('Eww180_2z.npy')

    #heights = [ 5, 30, 98 ];  
    heights = [ 5, 30 ];  
    numH = np.size(heights)
    ky = np.arange(0,ny/2)/2.0
    ky2 = np.arange(0,64)
    fig = plt.figure(figsize=(10,8))
    for j in range(0,numH):
        plt.subplot(1,numH,j+1)
        k = heights[j]
        plt.loglog(ky[1:], spec11[k,1:ny/2]*np.pi,'-',color='red',
                   label=r'$E_{uu}$')
        #plt.loglog(ky[1:], spec13[k+1,1:ny/2],'-', label=r'$E_{u'+vertvel+'}$')
        plt.loglog(ky[1:], spec22[k,1:ny/2]*np.pi,'-',color='blue', 
                   label=r'$E_{'+spanvel+spanvel+'}$')
        plt.loglog(ky[1:], spec33[k,1:ny/2]*np.pi,'-',color='green', 
                   label=r'$E_{'+vertvel+vertvel+'}$')
        if mkm:
            if j==0:
                plt.loglog(ky2, Euu180_0z,'o',color='red', label=r'$DNS, E_{uu}$')
                plt.loglog(ky2, Eww180_0z,'o',color='blue', 
                           label=r'$DNS, E_{'+spanvel+spanvel+'}$')
                plt.loglog(ky2, Evv180_0z,'o',color='green', 
                           label=r'$DNS, E_{'+vertvel+vertvel+'}$')
            if j==1:
                plt.loglog(ky2, Euu180_1z,'o',color='red', label=r'$DNS, E_{uu}$')
                plt.loglog(ky2, Eww180_1z,'o',color='blue', 
                           label=r'$DNS, E_{'+spanvel+spanvel+'}$')
                plt.loglog(ky2, Evv180_1z,'o',color='green',
                           label=r'$DNS, E_{'+vertvel+vertvel+'}$')
            if j==2:
                plt.loglog(ky2, Euu180_2z,'o',color='red', label=r'$DNS, E_{uu}$')
                plt.loglog(ky2, Eww180_2z,'o',color='blue', 
                           label=r'$DNS, E_{'+spanvel+spanvel+'}$')
                plt.loglog(ky2, Evv180_2z,'o',color='green',
                           label=r'$DNS, E_{'+vertvel+vertvel+'}$')

            
        plt.title(r'$ '+vert+'^{+} = $'+str(heights[j]));
        plt.xlabel(r'$ k_'+span+' $'); plt.tight_layout()
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
    lamY = Ly / ky
    
    kyE11 = ky * e11[:,0:ny/2]
    kyE13 = ky * e13[:,0:ny/2]
    kyE22 = ky * e22[:,0:ny/2]
    kyE33 = ky * e33[:,0:ny/2]

    if plot_wavelen:  # plot by wavelength
        xlab = '$ \lambda_'+span+' / \delta $';   ylab = '$'+vert+' / \delta $';
        myX, Z = np.meshgrid(lamY[1:], z[1:])
        tag = 'LAMY_'
    else:             # plot by wavenumber
        xlab = '$ \delta k_'+span+'$';   ylab = '$'+vert+' / \delta $';
        myX, Z = np.meshgrid(ky[1:], z[1:])
        tag = 'KY_'
    myFS = 18;
    numLevs = 100;

    fig = plt.figure(figsize=(12,4))

    ax = fig.add_subplot(2,2,1)
    levels = np.linspace(np.min(kyE11), np.max(kyE11), numLevs);
    kyE11plot = kyE11[1:,1:]
    cs = plt.contourf(myX, Z, kyE11plot, levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS); 
    plt.title(r'$k_'+span+r' E_{uu} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    np.max(kyE11plot)
    imax1,jmax1 = np.unravel_index(np.argmax(kyE11plot[:   ,0:10]),np.shape(kyE11plot[:,0:10]))  # split at ky=10
    imax2,jmax2 = np.unravel_index(np.argmax(kyE11plot[0:15,:]),np.shape(kyE11plot[0:15,:]))     # split at 15th z gridpoint
    print "peak1: ", kyE11plot[imax1,jmax1], myX[imax1,jmax1], Z[imax1,jmax1]
    print "peak2: ", kyE11plot[imax2,jmax2], myX[imax2,jmax2], Z[imax2,jmax2]

    arbNum = 50;
    plt.plot(np.ones(arbNum)*myX[imax1,jmax1], np.linspace(0,1,arbNum),'--',color='white')
    plt.plot(np.linspace(1,ny/2-1, arbNum),np.ones(arbNum)*Z[imax1,jmax1],'--',color='white')
    plt.plot(np.ones(arbNum)*myX[imax2,jmax2], np.linspace(0,1,arbNum),'--',color='white')
    plt.plot(np.linspace(1,ny/2-1, arbNum),np.ones(arbNum)*Z[imax2,jmax2],'--',color='white')
    #plt.plot(myX[imax1,jmax1], Z[imax1,jmax1],'o',color='white',markersize=7,markeredgewidth='2', fillstyle='none')
    #plt.plot(myX[imax2,jmax2+kysplit], Z[imax2,jmax2+kysplit],'o',color='white',markersize=7,markeredgewidth='2', fillstyle='none')
    plt.text(1.1,0.1,str(myX[imax1,jmax1])+', '+str(Z[imax1,jmax1]), fontsize=12, color='white')
    plt.text(10,0.1,str(myX[imax2,jmax2])+', '+str(Z[imax2,jmax2]), fontsize=12, color='white')
    
    ax = fig.add_subplot(2,2,2)
    levels = np.linspace(np.min(kyE13), np.max(kyE13), numLevs);
    cs = plt.contourf(myX, Z, kyE13[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS); 
    plt.title(r'$k_'+span+' E_{u'+vertvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()
    
    ax = fig.add_subplot(2,2,3)
    levels = np.linspace(np.min(kyE22), np.max(kyE22), numLevs);
    cs = plt.contourf(myX, Z, kyE22[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS); 
    plt.title(r'$k_'+span+' E_{'+spanvel+spanvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    imax1,jmax1 = np.unravel_index(np.argmax(kyE22),np.shape(kyE22))
    plt.plot(np.ones(arbNum)*myX[imax1,jmax1], np.linspace(0,1,arbNum),'--',color='white')
    plt.plot(np.linspace(1,ny/2-1, arbNum),np.ones(arbNum)*Z[imax1,jmax1],'--',color='white')
    #plt.plot(np.ones(arbNum)*28, np.linspace(0,1,arbNum),'--',color='white')
    #plt.plot(np.linspace(0,ny/2, arbNum),np.ones(arbNum)*0.04,'--',color='white')
    plt.text(1.1,0.02,str(myX[imax1,jmax1])+', '+str(Z[imax1,jmax1]), fontsize=12, color='white')

    ax = fig.add_subplot(2,2,4)
    levels = np.linspace(np.min(kyE33), np.max(kyE33), numLevs);
    cs = plt.contourf(myX, Z, kyE33[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS); 
    plt.title(r'$k_'+span+' E_{'+vertvel+vertvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    mySaveFig('sp1dky_'+tag, 0)

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
    lamX = Lx / kx
    
    kxE11 = kx * e11[:,0:nx/2]
    kxE13 = kx * e13[:,0:nx/2]
    kxE22 = kx * e22[:,0:nx/2]
    kxE33 = kx * e33[:,0:nx/2]

    if plot_wavelen:  # plot by wavelength
        xlab = '$ \lambda_x / \delta $';   ylab = '$'+vert+' / \delta $';
        myX, Z = np.meshgrid(lamX[1:], z[1:])
        tag = 'LAMX_'
    else:             # plot by wavenumber
        xlab = '$ \delta k_x$';   ylab = '$'+vert+' / \delta $';
        myX, Z = np.meshgrid(kx[1:], z[1:])
        tag = 'KX_'
    myFS = 18;
    numLevs = 100;

    fig = plt.figure(figsize=(12,4))

    ax = fig.add_subplot(2,2,1)
    levels = np.linspace(np.min(kxE11), np.max(kxE11), numLevs);
    kxE11plot = kxE11[1:,1:]
    cs = plt.contourf(myX, Z, kxE11plot, levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS); 
    plt.title(r'$k_x E_{uu} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    np.max(kxE11plot)
    imax1,jmax1 = np.unravel_index(np.argmax(kxE11plot[:   ,:]),np.shape(kxE11plot[:,:]))
    print "peak1: ", kxE11plot[imax1,jmax1], myX[imax1,jmax1], Z[imax1,jmax1]

    arbNum = 50;
    plt.plot(np.ones(arbNum)*myX[imax1,jmax1], np.linspace(0,1,arbNum),'--',color='white')
    plt.plot(np.linspace(1,nx/2-1, arbNum),np.ones(arbNum)*Z[imax1,jmax1],'--',color='white')
    plt.text(1.1,0.1,str(myX[imax1,jmax1])+', '+str(Z[imax1,jmax1]), fontsize=12, color='white')
    
    ax = fig.add_subplot(2,2,2)
    levels = np.linspace(np.min(kxE13), np.max(kxE13), numLevs);
    cs = plt.contourf(myX, Z, kxE13[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS); 
    plt.title(r'$k_x E_{u'+vertvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()
    
    ax = fig.add_subplot(2,2,3)
    levels = np.linspace(np.min(kxE22), np.max(kxE22), numLevs);
    cs = plt.contourf(myX, Z, kxE22[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS); 
    plt.title(r'$k_x E_{'+spanvel+spanvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    imax1,jmax1 = np.unravel_index(np.argmax(kxE22),np.shape(kxE22))
    plt.plot(np.ones(arbNum)*myX[imax1,jmax1], np.linspace(0,1,arbNum),'--',color='white')
    plt.plot(np.linspace(1,nx/2-1, arbNum),np.ones(arbNum)*Z[imax1,jmax1],'--',color='white')
    plt.text(1.1,0.08,str(myX[imax1,jmax1])+', '+str(Z[imax1,jmax1]), fontsize=12, color='white')

    ax = fig.add_subplot(2,2,4)
    levels = np.linspace(np.min(kxE33), np.max(kxE33), numLevs);
    cs = plt.contourf(myX, Z, kxE33[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS); 
    plt.title(r'$k_x E_{'+vertvel+vertvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    mySaveFig('sp1dkx_'+tag, 0)
    
if sp2d_plot_vert:
    sp2d_uu = np.load(datdir+'sp2d_uu.npy')
    sp2d_uw = np.load(datdir+'sp2d_uw.npy')
    sp2d_vv = np.load(datdir+'sp2d_vv.npy')
    sp2d_ww = np.load(datdir+'sp2d_ww.npy')

    suu = sp2d_uu[:, 0:ny/2, 0:nx/2]
    suu_sum = np.sum(suu, axis=2)  # sum over all kx

    kx = np.arange(0,nx/2)
    ky = np.arange(0,ny/2)
    lamY = Ly / ky
    lamX = Lx / kx
    
    LAMY, Z = np.meshgrid(lamY[1:], z[1:])
    KX, KY = np.meshgrid(kx[1:], ky[1:])

    fig = plt.figure(figsize=(12,8))
    xlab = '$ \lambda_'+span+' / \delta $';    ylab = '$'+vert+' / \delta $';    
    myFS = 18
    numLevs = 30;

    m1 = 0;     m2 = 1;     m3 = 2;     m4 = 3;

    ax = fig.add_subplot(2,2,1)
    #cs = plt.contourf(LAMY, Z, suu_sum[1:,1:])
    #levels = np.linspace(np.min(suu[:,:,m1]), np.max(suu[:,:,m1]), numLevs)
    cs = plt.contourf(LAMY, Z, suu[1:,1:,m1]) #, levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);  plt.ylabel(ylab, fontsize = myFS); 
    #plt.title(r'$E_{uu} all k_x$')
    plt.title( '$E_{uu}, k_x = '+str(m1)+'$' )
    cbar = plt.colorbar();   plt.tight_layout()

    ax = fig.add_subplot(2,2,2)
    levels = np.linspace(np.min(suu[:,:,m2]), np.max(suu[:,:,m2]), numLevs)
    cs = plt.contourf(LAMY, Z, suu[1:,1:,m2], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);  plt.ylabel(ylab, fontsize = myFS); 
    plt.title( '$E_{uu}, k_x = '+str(m2)+'$' )
    cbar = plt.colorbar();   plt.tight_layout()
     
    ax = fig.add_subplot(2,2,3)
    levels = np.linspace(np.min(suu[:,:,m3]), np.max(suu[:,:,m3]), numLevs)
    cs = plt.contourf(LAMY, Z, suu[1:,1:,m3], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);  plt.ylabel(ylab, fontsize = myFS); 
    plt.title( '$E_{uu}, k_x = '+str(m3)+'$' )
    cbar = plt.colorbar();   plt.tight_layout()
    
    ax = fig.add_subplot(2,2,4)
    levels = np.linspace(np.min(suu[:,:,m4]), np.max(suu[:,:,m4]), numLevs)
    cs = plt.contourf(LAMY, Z, suu[1:,1:,m4], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize = myFS);  plt.ylabel(ylab, fontsize = myFS); 
    plt.title( '$E_{uu}, k_x = '+str(m4)+'$' )
    cbar = plt.colorbar();   plt.tight_layout()

    mySaveFig('sp2d_vert', 0)

if sp2d_plot_horiz:
    sp2d_uu = np.load(datdir+'sp2d_uu.npy')
    sp2d_uw = np.load(datdir+'sp2d_uw.npy')
    sp2d_vv = np.load(datdir+'sp2d_vv.npy')
    sp2d_ww = np.load(datdir+'sp2d_ww.npy')

    xlab = '$ k_x / \delta $';    ylab = '$ k_'+span+' / \delta $';    
    myFS = 18;
    numLevs = 30;

    kxMax = nx/2;
    kyMax = ny/2;
    #kxMax =  5
    #kyMax =  20
    kx = np.arange(0,kxMax)
    ky = np.arange(0,kyMax)
    KX, KY = np.meshgrid(kx[1:], ky[1:])

    uu = sp2d_uu[:,1:kyMax,1:kxMax]
    vv = sp2d_vv[:,1:kyMax,1:kxMax]
    ww = sp2d_ww[:,1:kyMax,1:kxMax]
    uukx = uu * kx[1:]
    vvkx = vv * kx[1:]
    wwkx = ww * kx[1:]
    uukxky = uukx
    vvkxky = vvkx
    wwkxky = wwkx
    for j in range(1,np.size(ky)):
        uukxky[:,j-1,:] = uukxky[:,j-1,:] * ky[j]
        vvkxky[:,j-1,:] = vvkxky[:,j-1,:] * ky[j]
        wwkxky[:,j-1,:] = wwkxky[:,j-1,:] * ky[j]

    for k in range(1,nz,20):
        print 'k =', k
        sliceX = uukxky[k,:,:]
        sliceY = vvkxky[k,:,:]
        sliceZ = wwkxky[k,:,:]

        #sliceX =  sp2d_uu[k, 1:kyMax, 1:kxMax] * kx[1:]
        #sliceY =  sp2d_vv[k, 1:kyMax, 1:kxMax] * kx[1:]
        #sliceZ =  sp2d_ww[k, 1:kyMax, 1:kxMax] * kx[1:]
        #for j in range(1,np.size(ky)):
        #    sliceX[j-1,:] = sliceX[j-1,:] * ky[j]
        #    sliceY[j-1,:] = sliceY[j-1,:] * ky[j]
        #    sliceZ[j-1,:] = sliceZ[j-1,:] * ky[j]

        if localMax:
            minX = np.min(sliceX); maxX = np.max(sliceX)
            minY = np.min(sliceY); maxY = np.max(sliceY)
            minZ = np.min(sliceZ); maxZ = np.max(sliceZ)
            tag = '_loc_z'
        else:  # use global max instead
            minX = np.min(uukxky); maxX = np.max(uukxky)
            minY = np.min(vvkxky); maxY = np.max(vvkxky)
            minZ = np.min(wwkxky); maxZ = np.max(wwkxky)
            tag = '_glo_z'

        fig = plt.figure(figsize=(15,4))
        ax = fig.add_subplot(1,3,1)
        levels = np.linspace(minX, maxX, numLevs)
        cs = plt.contourf(KX, KY, sliceX, levels);
        ax.set_xscale('log');  ax.set_yscale('log')
        plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS);
        plt.title(r'$E_{uu}$'); cbar = plt.colorbar();    plt.tight_layout()

        ax = fig.add_subplot(1,3,2)
        levels = np.linspace(minY, maxY, numLevs)
        cs = plt.contourf(KX, KY, sliceY, levels);
        ax.set_xscale('log');  ax.set_yscale('log')
        plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS);
        plt.title(r'$E_{'+spanvel+spanvel+'}$'); 
        cbar = plt.colorbar();    plt.tight_layout()

        ax = fig.add_subplot(1,3,3)
        levels = np.linspace(minZ, maxZ, numLevs)
        cs = plt.contourf(KX, KY, sliceZ, levels);
        ax.set_xscale('log');  ax.set_yscale('log')
        plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS);
        plt.title(r'$E_{'+vertvel+vertvel+'}$'); 
        cbar = plt.colorbar();    plt.tight_layout()

        mySaveFig('sp2d' + tag + str(k)+'_', 0)
        plt.close()

    #ax = fig.add_subplot(2,2,2)
    #levels = np.linspace(np.min(slice2), np.max(slice2), numLevs)
    #cs = plt.contourf(KX, KY, slice2, levels);
    #ax.set_xscale('log');  ax.set_yscale('log')
    #plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS);
    #plt.title(r'$E_{uu}$')
    #cbar = plt.colorbar();    plt.tight_layout()
     
    #ax = fig.add_subplot(2,2,3)
    #levels = np.linspace(np.min(slice3), np.max(slice3), numLevs)
    #cs = plt.contourf(KX, KY, slice3, levels);
    #ax.set_xscale('log');  ax.set_yscale('log')
    #plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS);
    #plt.title(r'$E_{uu}$')
    #cbar = plt.colorbar();    plt.tight_layout()
    
    #ax = fig.add_subplot(2,2,4)
    #levels = np.linspace(np.min(slice4), np.max(slice4), numLevs)
    #cs = plt.contourf(KX, KY, slice4, levels);
    #ax.set_xscale('log');  ax.set_yscale('log')
    #plt.xlabel(xlab, fontsize = myFS);   plt.ylabel(ylab, fontsize = myFS);
    #plt.title(r'$E_{uu}$')
    #cbar = plt.colorbar();    plt.tight_layout()

    #mySaveFig('sp2d_horiz', 0)
    
    
if rs_plot:
    rs11Mean = np.load(datdir+'rs11Mean.npy')
    rs22Mean = np.load(datdir+'rs22Mean.npy')
    rs33Mean = np.load(datdir+'rs33Mean.npy')
    rs13Mean = np.load(datdir+'rs13Mean.npy')
    #rs23Mean = np.load(datdir+'rs23Mean.npy')
    #rs12Mean = np.load(datdir+'rs12Mean.npy')
    fig = plt.figure()
    #plt.subplot(2,1,1)
    L1=plt.plot(z, rs11Mean, 'o', color='b', 
                label = r'$[ u^{\prime} u^{\prime}]$')
    L2=plt.plot(z, rs22Mean, 'o', color='g', 
                label = r'$[ '+spanvel+'^{\prime} '+spanvel+'^{\prime}]$')
    L3=plt.plot(z, rs33Mean, 'o', color='r', 
                label = r'$['+vertvel+'^{\prime} '+vertvel+'^{\prime}]$')
    L4=plt.plot(z, rs13Mean, 'o', color='c', 
                label = r'$[ u^{\prime} '+vertvel+'^{\prime}]$')
    #L5=plt.plot(z, rs23Mean, 'o', color='m', 
    #            label = r'$[ '+spanvel+'^{\prime} '+vertvel+'^{\prime}]$')
    #L6=plt.plot(z, rs12Mean, 'o', color='k', 
    #            label = r'$[ u^{\prime} '+spanvel+'^{\prime}]$')
    plt.legend(loc='upper right',ncol=2)
    plt.xlabel(r'$'+vert+' / H $', fontsize=18); 
    plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.ylim([-2, 20])
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
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
    L2=plt.plot(z, rs22Mean, 'o', color='g', 
                label = r'$[ '+spanvel+'^{\prime} '+spanvel+'^{\prime}]$')
    L3=plt.plot(z, rs33Mean, 'o', color='r', 
                label = r'$[ '+vertvel+'^{\prime} '+vertvel+'^{\prime}]$')
    L4=plt.plot(z, -1.0*rs13Mean, 'o', color='c', 
                label = r'$[ -u^{\prime} '+vertvel+'^{\prime}]$')
    #L5=plt.plot(z, rs23Mean, 'o', color='m', 
    #            label = r'$[ '+spanvel+'^{\prime} '+vertvel+'^{\prime}]$')
    #L6=plt.plot(z, rs12Mean, 'o', color='k', 
    #            label = r'$[ u^{\prime} '+spanvel+'^{\prime}]$')
    plt.legend(loc='upper right',ncol=2)
    plt.xlabel(r'$ '+vert+' / H $', fontsize=18); 
    plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.tight_layout()
    mySaveFig('rs2_', 0)

    if vel2_rs_plot:
            fig = plt.figure()
            vvMean = np.load(datdir+'vvMean.npy')
            wwMean = np.load(datdir+'wwMean.npy')
            uwMean = np.load(datdir+'uwMean.npy')
            vwMean = np.load(datdir+'vwMean.npy')
            uvMean = np.load(datdir+'uvMean.npy')
            #L1=plt.plot(z, rs11Mean, 'o', color='b', 
            #            label = r'$[ u^{\prime} u^{\prime}]$')
            L2=plt.plot(z, rs22Mean, 'o', color='g', 
                        label = r'$[ '+spanvel+'^{\prime} '+spanvel+'^{\prime}]$')
            L3=plt.plot(z, rs33Mean, 'o', color='r', 
                        label = r'$[ '+vertvel+'^{\prime} '+vertvel+'^{\prime}]$')
            L4=plt.plot(z, -1*rs13Mean, 'o', color='c', 
                        label = r'$[ u^{\prime} '+vertvel+'^{\prime}]$')
            L5=plt.plot(z, rs23Mean, 'o', color='m', 
                        label = r'$[ '+spanvel+'^{\prime} '+vertvel+'^{\prime}]$')
            L6=plt.plot(z, rs12Mean, 'o', color='k', 
                        label = r'$[ u^{\prime} '+spanvel+'^{\prime}]$')
            #L7=plt.plot(z, uuMean, color='b', label = r'$[ uu ]$')
            L8=plt.plot(z, vvMean, linewidth=3, color='g', 
                        label = r'$[ '+spanvel+spanvel+' ]$')
            L9=plt.plot(z, wwMean, linewidth=3, color='r', 
                        label = r'$[ '+vertvel+vertvel+' ]$')
            L10=plt.plot(z, -1*uwMean, linewidth=3, color='c', 
                         label = r'$[ u'+vertvel+' ]$')
            L11=plt.plot(z, vwMean, linewidth=3, color='m', 
                         label = r'$[ '+spanvel+vertvel+' ]$')
            L12=plt.plot(z, uvMean, linewidth=3, color='k', 
                         label = r'$[ u'+spanvel+' ]$')
            plt.legend(loc='lower right',ncol=2)
            plt.xlabel(r'$ '+vert+' / H $', fontsize=18); 
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
            plt.plot(z,vvMean,'s',label=spanvel+spanvel+'Mean')
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
            plt.plot(z,wwMean,'s',label=vertvel+vertvel+'Mean')
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
        plt.xlabel(r'$ x / H $', fontsize=18); 
        plt.ylabel(r'$'+span+' / H $', fontsize=18); 
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
        plt.xlabel(r'$'+span+' / H $', fontsize=18); 
        plt.ylabel(r'$'+vert+' / H $', fontsize=18); 
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
