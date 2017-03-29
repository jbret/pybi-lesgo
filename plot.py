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

matplotlib.rcParams['lines.linewidth'] = 1.5
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['lines.markersize'] = 8
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['xtick.major.pad'] = 10
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['ytick.major.pad'] = 10
matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['axes.labelsize'] = 34 # 34 for paper
matplotlib.rcParams['legend.fontsize'] = 22
matplotlib.rcParams.update({'figure.autolayout': True})

# For font types (Journal does not accept type 3 fonts)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

RNL_branch = 1;    devel_branch = 0;

vel_avg_plot    = 1;
uXMean_plot     = 0;
tau_plot        = 1;
nu_t_plot       = 0;
spanSpec_plot   = 0;
sp1dky_plot     = 1;  plot_wavelen = 0;  # by wavelength or wavenumber
sp1dkx_plot     = 0;
spvort_plot     = 0;  vort_components = 1;   vort_more = 0;  integrate_spectra = 0; 
sp2d_plot_vert  = 0;  # WARNING: must test plot labels in LES case for sp2d plots (both vert and horiz)
sp2d_plot_horiz = 0;  localMax = 1  # if 0 then uses global max
rs_plot         = 1;
vel2_rs_plot    = 0;
snap_plot_xy    = 0;   thisSnap = 250300;  # on uv-grid
snap_plot_yz    = 0;
energy_bal      = 0;
test_plot       = 0;
#snap_plot       = 1;  thisSnap = 250300;  # on uv-grid

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
    plt.semilogx(z, 1/0.4*np.log(z/.0000125), '-k', 
                 label=r'$1/\kappa \ \mathrm{log}('+vert+'/'+vert+'_{0})$')
    #plt.semilogx(z/(1./180), 1/0.41*np.log(z/(1./180))+5.0, '-k', label=r'$1/\kappa \ \mathrm{log}(z^{+})+B$')
    plt.xlabel('$'+vert+' / H $');
    #plt.xlabel('$'+vert+'^{+} $');
    plt.ylabel('$[\ u / u_*]$'); 
    #plt.xlim([.02 ,1.1])
    #plt.xlim([1, 1000])
    plt.text(.2,6,r'$ \kappa = 0.4,\ '+vert+'_{0} = 10^{-4} $', fontsize=16)
    #plt.text(.06,1,r'$ \kappa = 0.4,\ '+vert+'_{0} = 10^{-4} $', fontsize=14)
    plt.xticks([10**(-2), 10**(-1), 10**0])
    plt.yticks([0,5,10,15,20,25])

    if mkm:
        yp180 = np.load('yp180.npy')
        Umean180 = np.load('Umean180.npy')
        plt.semilogx(yp180, Umean180,'o',label=r'$MKM99, DNS, Re_{\tau} = 180 $')

    plt.text(110,2.5,r'$ \kappa = 0.41,\ B = 5.0 $', fontsize=14)
    plt.legend(loc='lower right')
    plt.tight_layout()

    mySaveFig('mvp_', 0)

if uXMean_plot:
    scale = 3.0;
    fig = plt.figure(figsize=(scale*Ly,scale*Lz))
    Y, Z = np.meshgrid(y, z)
    uXMean = np.load(datdir+'uXMean.npy')
    cs = plt.contourf(Y, Z, uXMean[:,:], vmin=0, vmax=17)
    cbar = plt.colorbar()
    plt.xlabel('$'+span+' / H $'); 
    plt.ylabel('$'+vert+' / H $');
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
    plt.xlabel('$ \mathrm{Stress} $'); 
    plt.ylabel('$'+vert+' / H $')
    plt.legend()
    plt.tight_layout()
    mySaveFig('tau_', 0)

if nu_t_plot:
    nu_tMean = np.load(datdir+'nu_tMean.npy')
    fig = plt.figure(figsize=(14,7))
    pf = 1.0/180.0    #450, 1050
    plt.subplot(1,2,1)
    plt.plot(nu_tMean, z, '-o', color='g', label = r'$[ \nu_T ]$')
    plt.plot(pf*z**(-1/3.0), z, '-', color='r', label = r'$c*z^{-1/3}$')
    # note below the 'r' needed before label (in order to render the \nu)
    plt.xlabel(r'$[ \nu_T ]$'); 
    plt.ylabel('$'+vert+' / H $')
    plt.legend()
    plt.tight_layout()
    plt.subplot(1,2,2)
    plt.loglog(nu_tMean, z, '-o', color='g', label = r'$[ \nu_T ]$')
    plt.loglog(pf*z**(-1/3.0), z, '-', color='r', label = r'$c*z^{-1/3}$')
    # note below the 'r' needed before label (in order to render the \nu)
    plt.xlabel(r'$[ \nu_T ]$');
    plt.ylabel('$'+vert+' / H $')
    plt.legend()
    plt.tight_layout()
    mySaveFig('nu_t_', 0)

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
    subplot_fs  = 18;
    subplot_fs2 = 14;
    numLevs = 100;

    fig = plt.figure(figsize=(12,4))

    ax = fig.add_subplot(2,2,1)
    levels = np.linspace(np.min(kyE11), np.max(kyE11), numLevs);
    kyE11plot = kyE11[1:,1:]
    cs = plt.contourf(myX, Z, kyE11plot, levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize=subplot_fs);   plt.ylabel(ylab, fontsize=subplot_fs);
    plt.xticks(fontsize=subplot_fs2);  plt.yticks([10**(-1), 10**0],fontsize=subplot_fs2) 
    plt.title(r'$k_'+span+r' E_{uu} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    np.max(kyE11plot)
    imax1,jmax1 = np.unravel_index(np.argmax(kyE11plot[:   ,0:10]),np.shape(kyE11plot[:,0:10]))  # split at ky=10
    imax2,jmax2 = np.unravel_index(np.argmax(kyE11plot[0:15,:]),np.shape(kyE11plot[0:15,:]))     # split at 15th z gridpoint
    print "peak1: ", kyE11plot[imax1,jmax1], myX[imax1,jmax1], Z[imax1,jmax1]
    print "peak2: ", kyE11plot[imax2,jmax2], myX[imax2,jmax2], Z[imax2,jmax2]

    #arbNum = 50;
    #plt.plot(np.ones(arbNum)*myX[imax1,jmax1], np.linspace(0.01,1,arbNum),'--',color='white')
    #plt.plot(np.linspace(1,ny/2-1, arbNum),np.ones(arbNum)*Z[imax1,jmax1],'--',color='white')
    #plt.plot(np.ones(arbNum)*myX[imax2,jmax2], np.linspace(0.01,1,arbNum),'--',color='white')
    #plt.plot(np.linspace(1,ny/2-1, arbNum),np.ones(arbNum)*Z[imax2,jmax2],'--',color='white')
    ##plt.plot(myX[imax1,jmax1], Z[imax1,jmax1],'o',color='white',markersize=7,markeredgewidth='2', fillstyle='none')
    ##plt.plot(myX[imax2,jmax2+kysplit], Z[imax2,jmax2+kysplit],'o',color='white',markersize=7,markeredgewidth='2', fillstyle='none')
    #plt.text(1.1,0.1,str(myX[imax1,jmax1])+', '+str(Z[imax1,jmax1]), fontsize=12, color='white')
    #plt.text(10,0.1,str(myX[imax2,jmax2])+', '+str(Z[imax2,jmax2]), fontsize=12, color='white')
    
    ax = fig.add_subplot(2,2,2)
    levels = np.linspace(np.min(kyE13), np.max(kyE13), numLevs);
    cs = plt.contourf(myX, Z, kyE13[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize=subplot_fs);   plt.ylabel(ylab, fontsize=subplot_fs); 
    plt.xticks(fontsize=subplot_fs2);  plt.yticks([10**(-1), 10**0],fontsize=subplot_fs2) 
    plt.title(r'$k_'+span+' E_{u'+vertvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()
    
    ax = fig.add_subplot(2,2,3)
    levels = np.linspace(np.min(kyE22), np.max(kyE22), numLevs);
    cs = plt.contourf(myX, Z, kyE22[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize=subplot_fs);   plt.ylabel(ylab, fontsize=subplot_fs); 
    plt.xticks(fontsize=subplot_fs2);  plt.yticks([10**(-1), 10**0],fontsize=subplot_fs2) 
    plt.title(r'$k_'+span+' E_{'+spanvel+spanvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    #imax1,jmax1 = np.unravel_index(np.argmax(kyE22),np.shape(kyE22))
    #plt.plot(np.ones(arbNum)*myX[imax1,jmax1], np.linspace(0.01,1,arbNum),'--',color='white')
    #plt.plot(np.linspace(1,ny/2-1, arbNum),np.ones(arbNum)*Z[imax1,jmax1],'--',color='white')
    ##plt.plot(np.ones(arbNum)*28, np.linspace(0.01,1,arbNum),'--',color='white')
    ##plt.plot(np.linspace(0,ny/2, arbNum),np.ones(arbNum)*0.04,'--',color='white')
    #plt.text(1.1,0.02,str(myX[imax1,jmax1])+', '+str(Z[imax1,jmax1]), fontsize=12, color='white')

    ax = fig.add_subplot(2,2,4)
    levels = np.linspace(np.min(kyE33), np.max(kyE33), numLevs);
    cs = plt.contourf(myX, Z, kyE33[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab, fontsize=subplot_fs);   plt.ylabel(ylab, fontsize=subplot_fs); 
    plt.xticks(fontsize=subplot_fs2);  plt.yticks([10**(-1), 10**0],fontsize=subplot_fs2) 
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
    numLevs = 100;

    fig = plt.figure(figsize=(12,4))

    ax = fig.add_subplot(2,2,1)
    levels = np.linspace(np.min(kxE11), np.max(kxE11), numLevs);
    kxE11plot = kxE11[1:,1:]
    cs = plt.contourf(myX, Z, kxE11plot, levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab);   plt.ylabel(ylab); 
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
    plt.xlabel(xlab);   plt.ylabel(ylab); 
    plt.title(r'$k_x E_{u'+vertvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()
    
    ax = fig.add_subplot(2,2,3)
    levels = np.linspace(np.min(kxE22), np.max(kxE22), numLevs);
    cs = plt.contourf(myX, Z, kxE22[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab);   plt.ylabel(ylab); 
    plt.title(r'$k_x E_{'+spanvel+spanvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    #imax1,jmax1 = np.unravel_index(np.argmax(kxE22),np.shape(kxE22))
    #plt.plot(np.ones(arbNum)*myX[imax1,jmax1], np.linspace(0,1,arbNum),'--',color='white')
    #plt.plot(np.linspace(1,nx/2-1, arbNum),np.ones(arbNum)*Z[imax1,jmax1],'--',color='white')
    #plt.text(1.1,0.08,str(myX[imax1,jmax1])+', '+str(Z[imax1,jmax1]), fontsize=12, color='white')

    ax = fig.add_subplot(2,2,4)
    levels = np.linspace(np.min(kxE33), np.max(kxE33), numLevs);
    cs = plt.contourf(myX, Z, kxE33[1:,1:], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab);   plt.ylabel(ylab); 
    plt.title(r'$k_x E_{'+vertvel+vertvel+r'} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()

    mySaveFig('sp1dkx_'+tag, 0)

if spvort_plot:
    # NOTE: Figures created through the pyplot interface (`matplotlib.pyplot.figure`) 
    # are retained until explicitly closed and may consume too much memory.
    # Could have problems if over 20 figures are opened.
    # --------------------------------------------------------
    # square of vorticity component recorded in physical space
    spvort_vortx = np.load(datdir+'spvort_vortx.npy')
    spvort_vorty = np.load(datdir+'spvort_vorty.npy')
    spvort_vortz = np.load(datdir+'spvort_vortz.npy')
    # sum of squared vorticity components recorded in physical space
    # i.e., the square of the vorticity magnitude
    spvort_vortp = np.load(datdir+'spvort_vortp.npy')
    # fourier component multiplied by its complex conj (kx space)
    spvort_vortsx = np.load(datdir+'spvort_vortsx.npy')
    spvort_vortsy = np.load(datdir+'spvort_vortsy.npy')
    spvort_vortsz = np.load(datdir+'spvort_vortsz.npy')
    # fourier components (of vorticity magnitude) mult by its complex conj (kx space)
    # --> square root was taken in calculating vorticity magnitude
    spvort_vorts = np.load(datdir+'spvort_vorts.npy')

    # spanwise average
    ep = np.mean(spvort_vortp[:,:,:], axis=1)
    ep_cross = np.mean(spvort_vortp[:,:,:], axis=2)
    ep_cross2 = np.mean(spvort_vortp[:,:,:], axis=0)
    ex = np.mean(spvort_vortx[:,:,:], axis=1)
    ex_cross = np.mean(spvort_vortx[:,:,:], axis=2)
    ex_cross2 = np.mean(spvort_vortx[:,:,:], axis=0)
    ey = np.mean(spvort_vorty[:,:,:], axis=1)
    ey_cross = np.mean(spvort_vorty[:,:,:], axis=2)
    ey_cross2 = np.mean(spvort_vorty[:,:,:], axis=0)
    ez = np.mean(spvort_vortz[:,:,:], axis=1)
    ez_cross = np.mean(spvort_vortz[:,:,:], axis=2)
    ez_cross2 = np.mean(spvort_vortz[:,:,:], axis=0)
    
    es = np.mean(spvort_vorts[:,:,:], axis=1)
    esx = np.mean(spvort_vortsx[:,:,:], axis=1)
    esy = np.mean(spvort_vortsy[:,:,:], axis=1)
    esz = np.mean(spvort_vortsz[:,:,:], axis=1)

    kxMax = nx/2
    #kxMax = 15
    kx = np.arange(0,kxMax)
    lamX = Lx / kx
    
    kxEs  = kx * es [:, 0:kxMax]
    kxEsx = kx * esx[:, 0:kxMax]
    kxEsy = kx * esy[:, 0:kxMax]
    kxEsz = kx * esz[:, 0:kxMax]
    #kxEsSum = np.zeros(np.shape(kxEs))
    #kxEsSum = kxEsx + kxEsy + kxEsz
    kxEsSum = esx[:,0:kxMax] + esy[:,0:kxMax] + esz[:,0:kxMax]
    #sumHeight = np.sum(kxEsSum,1)
    #tot = np.sum(sumHeight)
    #kxEsSum = kxEsSum / tot
    kxEsSum = kx * kxEsSum

    kxCap = kxMax

    # cannot exactly recover the physical space version because
    # the square root was taken when calculating spvort_vorts
    ints_ = spvort_vorts[:,:,0:kxCap]
    ints_[:,:,1:] = ints_[:,:,1:] * 2.0
    ints_ = np.sum(ints_[:,:,:], axis=2) # integrate over kx
    ints  = np.mean(ints_[:,:], axis=1)  # average over span

    intx_ = spvort_vortsx[:,:,0:kxCap]
    intx_[:,:,1:] = intx_[:,:,1:] * 2.0
    intx_ = np.sum(intx_[:,:,:], axis=2) # integrate over kx
    intx = np.mean(intx_[:,:], axis=1)   # average over span

    inty_ = spvort_vortsy[:,:,0:kxCap]
    inty_[:,:,1:] = inty_[:,:,1:] * 2.0
    inty_ = np.sum(inty_[:,:,:], axis=2) # integrate over kx
    inty = np.mean(inty_[:,:], axis=1)   # average over span

    intz_ = spvort_vortsz[:,:,0:kxCap]
    intz_[:,:,1:] = intz_[:,:,1:] * 2.0
    intz_ = np.sum(intz_[:,:,:], axis=2) # integrate over kx
    intz = np.mean(intz_[:,:], axis=1)   # average over span

    intTot = intx + inty + intz
    
    for jz in range(1,nz):
        kxEsSum[jz,:] = kxEsSum[jz,:]  #/ intTot[jz]
        kxEsx[jz,:] = kxEsx[jz,:]      #/ intx[jz]
        kxEsy[jz,:] = kxEsy[jz,:]      #/ inty[jz]
        kxEsz[jz,:] = kxEsz[jz,:]      #/ intz[jz]

    A = 0.86;  B = 1.0;  # for re-scaling as in Jimenez review paper
    kxEsNorm = np.zeros(np.shape(kxEs))
    f = kxEs
    #f = kxEsSum
    for jz in range(1,nz):
        myMin = np.min( f[jz,1:] )
        myMax = np.max( f[jz,1:] )
        # rescale [myMin, myMax] into [A, B]
        kxEsNorm[jz,:] = A + ( f[jz,:]-myMin )*(B-A)/(myMax-myMin)

    if plot_wavelen:  # plot by wavelength
        xlab = '$ \lambda_x / \delta $';   ylab = '$'+vert+' / \delta $';
        myX, Z = np.meshgrid(lamX[1:], z[1:])
        tag = 'LAMX_'+'kx'+str(kxMax)+'_'
    else:             # plot by wavenumber
        xlab = '$ \delta k_x$';   ylab = '$'+vert+' / \delta $';
        myX, Z = np.meshgrid(kx[1:], z[1:])
        tag = 'KX_'+'kx'+str(kxMax)+'_'
    numLevs = 100;

    fig = plt.figure()
    kxEsplot = kxEsSum[1:,1:]
    levels = np.linspace(np.min(kxEsplot), np.max(kxEsplot), numLevs);
    cs = plt.contourf(myX, Z, kxEsplot, levels);
    plt.yscale('log'); plt.xscale('log')
    plt.xlabel(xlab);   plt.ylabel(ylab); 
    #plt.title(r'$k_x E_{uu} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()
    mySaveFig('vorts_'+tag, 0); plt.close()
    #np.save('myX32',myX)
    #np.save('Z32',Z)
    #np.save('kxEsplot32',kxEsplot)

    fig = plt.figure()
    kxEsNormplot = kxEsNorm[1:,1:]
    levels = np.linspace(np.min(kxEsNormplot), np.max(kxEsNormplot), numLevs);
    cs = plt.contourf(myX, Z, kxEsNormplot, levels);
    plt.yscale('log'); plt.xscale('log')
    plt.xlabel(xlab);   plt.ylabel(ylab); 
    #plt.title(r'$k_x E_{uu} / u_{\tau}^2 $')
    cbar = plt.colorbar();    plt.tight_layout()
    mySaveFig('vortsNormMagOnly_'+tag, 0); plt.close()


    if vort_components == 1:
        fig = plt.figure()
        kxEsxplot = kxEsx[1:,1:]
        levels = np.linspace(np.min(kxEsxplot), np.max(kxEsxplot), numLevs);
        cs = plt.contourf(myX, Z, kxEsxplot, levels);
        plt.yscale('log'); plt.xscale('log')
        plt.xlabel(xlab);   plt.ylabel(ylab); 
        #plt.title(r'$k_x E_{uu} / u_{\tau}^2 $')
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortsx_'+tag, 0); plt.close()
        #np.save('kxEsxplot32',kxEsxplot)

        fig = plt.figure()
        kxEsyplot = kxEsy[1:,1:]
        levels = np.linspace(np.min(kxEsyplot), np.max(kxEsyplot), numLevs);
        cs = plt.contourf(myX, Z, kxEsyplot, levels);
        plt.yscale('log'); plt.xscale('log')
        plt.xlabel(xlab);   plt.ylabel(ylab); 
        #plt.title(r'$k_x E_{uu} / u_{\tau}^2 $')
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortsy_'+tag, 0); plt.close()
        #np.save('kxEsyplot32',kxEsyplot)        

        fig = plt.figure()
        kxEszplot = kxEsz[1:,1:]
        levels = np.linspace(np.min(kxEszplot), np.max(kxEszplot), numLevs);
        cs = plt.contourf(myX, Z, kxEszplot, levels);
        imax1,jmax1=np.unravel_index(np.argmax(kxEszplot[:,:]),np.shape(kxEszplot[:,:]))
        print "peak: ", kxEszplot[imax1,jmax1], myX[imax1,jmax1], kxEszplot[imax1,jmax1]
        arbNum = 50;
        plt.plot(np.ones(arbNum)*myX[imax1,jmax1], np.linspace(z[1],z[-1],arbNum),'--',color='white')
        plt.plot(np.linspace(1,kxMax-1, arbNum),np.ones(arbNum)*Z[imax1,jmax1],'--',color='white')
        plt.text(1.1,0.8,r'$k_x='+str(myX[imax1,jmax1])+'$', fontsize=12, color='white')
        plt.yscale('log'); plt.xscale('log')
        plt.xlabel(xlab);   plt.ylabel(ylab); 
        #plt.title(r'$k_x E_{uu} / u_{\tau}^2 $')
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortsz_'+tag, 0); plt.close()
        #np.save('kxEszplot32',kxEszplot)

    if vort_more == 1:
        fig = plt.figure()
        X, Z = np.meshgrid(x, z)
        cs = plt.contourf(X[:,:], Z[:,:], ep[:,:]);
        plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortp_'+tag, 0); plt.close()

        fig = plt.figure()
        Y, Z = np.meshgrid(y, z)
        cs = plt.contourf(Y[:,:], Z[:,:], ep_cross[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortp_cross_'+tag, 0); plt.close()

        fig = plt.figure()
        X, Y = np.meshgrid(x, y)
        cs = plt.contourf(X[:,:], Y[:,:], ep_cross2[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortp_cross2_'+tag, 0); plt.close()

        fig = plt.figure()
        X, Z = np.meshgrid(x, z)
        cs = plt.contourf(X[:,:], Z[:,:], ex[:,:]);
        plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortx_'+tag, 0); plt.close()

        fig = plt.figure()
        Y, Z = np.meshgrid(y, z)
        cs = plt.contourf(Y[:,:], Z[:,:], ex_cross[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortx_cross_'+tag, 0); plt.close()

        fig = plt.figure()
        X, Y = np.meshgrid(x, y)
        cs = plt.contourf(X[:,:], Y[:,:], ex_cross2[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortx_cross2_'+tag, 0); plt.close()

        fig = plt.figure()
        X, Z = np.meshgrid(x, z)
        cs = plt.contourf(X[:,:], Z[:,:], ey[:,:]);
        plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vorty_'+tag, 0); plt.close()

        fig = plt.figure()
        Y, Z = np.meshgrid(y, z)
        cs = plt.contourf(Y[:,:], Z[:,:], ey_cross[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vorty_cross_'+tag, 0); plt.close()

        fig = plt.figure()
        X, Y = np.meshgrid(x, y)
        cs = plt.contourf(X[:,:], Y[:,:], ey_cross2[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vorty_cross2_'+tag, 0); plt.close()

        fig = plt.figure()
        X, Z = np.meshgrid(x, z)
        cs = plt.contourf(X[:,:], Z[:,:], ez[:,:]);
        plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortz_'+tag, 0); plt.close()

        fig = plt.figure()
        Y, Z = np.meshgrid(y, z)
        cs = plt.contourf(Y[:,:], Z[:,:], ez_cross[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortz_cross_'+tag, 0); plt.close()

        fig = plt.figure()
        X, Y = np.meshgrid(x, y)
        cs = plt.contourf(X[:,:], Y[:,:], ez_cross2[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('vortz_cross2_'+tag, 0); plt.close()
    if integrate_spectra:
        # Sum spectra and compare
        #compx_ = esz[:,0:nx/2]
        #compx_[:,1:nx/2] = compx_[:,1:nx/2] + compx_[:,1:nx/2]
        #compx = np.sum(compx_, axis=1)
        #fig = plt.figure()
        #comp  = np.mean(ez,axis=1)
        #plt.plot(z, comp,'sg')
        #plt.plot(z, compx,'or')
        
        fig = plt.figure()
        Y, Z = np.meshgrid(y, z)
        ax = fig.add_subplot(1,2,1)
        cs = plt.contourf(Y[:,:], Z[:,:], ep_cross[:,:]);
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        ax = fig.add_subplot(1,2,2)
        cs = plt.contourf( Y[:,:], Z[:,:], ints_[:,:] );
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('comp_',0);     plt.close()
        ## ^^ won't match exactly because of square root
        ## but the individual components (below) will match

        fig = plt.figure()
        Y, Z = np.meshgrid(y, z)
        ax = fig.add_subplot(1,2,1)
        cs = plt.contourf( Y[:,:], Z[:,:], ex_cross[:,:] );
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        ax = fig.add_subplot(1,2,2)
        cs = plt.contourf( Y[:,:], Z[:,:], intx_[:,:] );
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('compx_',0);    plt.close()
        print "difference: ", np.sum( ex_cross[:,:] - intx_[:,:] )

        fig = plt.figure()
        Y, Z = np.meshgrid(y, z)
        ax = fig.add_subplot(1,2,1)
        cs = plt.contourf( Y[:,:], Z[:,:], ey_cross[:,:] );
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        ax = fig.add_subplot(1,2,2)
        cs = plt.contourf( Y[:,:], Z[:,:], inty_[:,:] );
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('compy_',0);    plt.close()
        print "difference: ", np.sum( ey_cross[:,:] - inty_[:,:] )

        fig = plt.figure()
        Y, Z = np.meshgrid(y, z)
        ax = fig.add_subplot(1,2,1)
        cs = plt.contourf( Y[:,:], Z[:,:], ez_cross[:,:] );
        #plt.xlabel(xlab);   plt.ylabel(ylab); 
        cbar = plt.colorbar();    plt.tight_layout()
        ax = fig.add_subplot(1,2,2)
        cs = plt.contourf( Y[:,:], Z[:,:], intz_[:,:] );
        cbar = plt.colorbar();    plt.tight_layout()
        mySaveFig('compz_',0);    plt.close()
        print "difference: ", np.sum( ez_cross[:,:] - intz_[:,:] )

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
    plt.xlabel(xlab);  plt.ylabel(ylab); 
    #plt.title(r'$E_{uu} all k_x$')
    plt.title( '$E_{uu}, k_x = '+str(m1)+'$' )
    cbar = plt.colorbar();   plt.tight_layout()

    ax = fig.add_subplot(2,2,2)
    levels = np.linspace(np.min(suu[:,:,m2]), np.max(suu[:,:,m2]), numLevs)
    cs = plt.contourf(LAMY, Z, suu[1:,1:,m2], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab);  plt.ylabel(ylab); 
    plt.title( '$E_{uu}, k_x = '+str(m2)+'$' )
    cbar = plt.colorbar();   plt.tight_layout()
     
    ax = fig.add_subplot(2,2,3)
    levels = np.linspace(np.min(suu[:,:,m3]), np.max(suu[:,:,m3]), numLevs)
    cs = plt.contourf(LAMY, Z, suu[1:,1:,m3], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab);  plt.ylabel(ylab); 
    plt.title( '$E_{uu}, k_x = '+str(m3)+'$' )
    cbar = plt.colorbar();   plt.tight_layout()
    
    ax = fig.add_subplot(2,2,4)
    levels = np.linspace(np.min(suu[:,:,m4]), np.max(suu[:,:,m4]), numLevs)
    cs = plt.contourf(LAMY, Z, suu[1:,1:,m4], levels);
    ax.set_xscale('log');  ax.set_yscale('log')
    plt.xlabel(xlab);  plt.ylabel(ylab); 
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
        plt.xlabel(xlab);   plt.ylabel(ylab);
        plt.title(r'$E_{uu}$'); cbar = plt.colorbar();    plt.tight_layout()

        ax = fig.add_subplot(1,3,2)
        levels = np.linspace(minY, maxY, numLevs)
        cs = plt.contourf(KX, KY, sliceY, levels);
        ax.set_xscale('log');  ax.set_yscale('log')
        plt.xlabel(xlab);   plt.ylabel(ylab);
        plt.title(r'$E_{'+spanvel+spanvel+'}$'); 
        cbar = plt.colorbar();    plt.tight_layout()

        ax = fig.add_subplot(1,3,3)
        levels = np.linspace(minZ, maxZ, numLevs)
        cs = plt.contourf(KX, KY, sliceZ, levels);
        ax.set_xscale('log');  ax.set_yscale('log')
        plt.xlabel(xlab);   plt.ylabel(ylab);
        plt.title(r'$E_{'+vertvel+vertvel+'}$'); 
        cbar = plt.colorbar();    plt.tight_layout()

        mySaveFig('sp2d' + tag + str(k)+'_', 0)
        plt.close()

    #ax = fig.add_subplot(2,2,2)
    #levels = np.linspace(np.min(slice2), np.max(slice2), numLevs)
    #cs = plt.contourf(KX, KY, slice2, levels);
    #ax.set_xscale('log');  ax.set_yscale('log')
    #plt.xlabel(xlab);   plt.ylabel(ylab);
    #plt.title(r'$E_{uu}$')
    #cbar = plt.colorbar();    plt.tight_layout()
     
    #ax = fig.add_subplot(2,2,3)
    #levels = np.linspace(np.min(slice3), np.max(slice3), numLevs)
    #cs = plt.contourf(KX, KY, slice3, levels);
    #ax.set_xscale('log');  ax.set_yscale('log')
    #plt.xlabel(xlab);   plt.ylabel(ylab);
    #plt.title(r'$E_{uu}$')
    #cbar = plt.colorbar();    plt.tight_layout()
    
    #ax = fig.add_subplot(2,2,4)
    #levels = np.linspace(np.min(slice4), np.max(slice4), numLevs)
    #cs = plt.contourf(KX, KY, slice4, levels);
    #ax.set_xscale('log');  ax.set_yscale('log')
    #plt.xlabel(xlab);   plt.ylabel(ylab);
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
    plt.xlabel(r'$'+vert+' / H $'); 
    plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $')
    plt.ylim([-2, 20])
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
    #plt.xlabel(r'$ z^{+} $'); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $')
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
    plt.legend(loc='upper right',ncol=1)
    plt.xlabel(r'$ '+vert+' / H $'); 
    plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $')
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
            #L5=plt.plot(z, rs23Mean, 'o', color='m', 
            #            label = r'$[ '+spanvel+'^{\prime} '+vertvel+'^{\prime}]$')
            #L6=plt.plot(z, rs12Mean, 'o', color='k', 
            #            label = r'$[ u^{\prime} '+spanvel+'^{\prime}]$')
            #L7=plt.plot(z, uuMean, color='b', label = r'$[ uu ]$')
            L8=plt.plot(z, vvMean, linewidth=3, color='g', 
                        label = r'$[ '+spanvel+spanvel+' ]$')
            L9=plt.plot(z, wwMean, linewidth=3, color='r', 
                        label = r'$[ '+vertvel+vertvel+' ]$')
            L10=plt.plot(z, -1*uwMean, linewidth=3, color='c', 
                         label = r'$[ u'+vertvel+' ]$')
            #L11=plt.plot(z, vwMean, linewidth=3, color='m', 
            #             label = r'$[ '+spanvel+vertvel+' ]$')
            #L12=plt.plot(z, uvMean, linewidth=3, color='k', 
            #             label = r'$[ u'+spanvel+' ]$')
            plt.legend(loc='lower right',ncol=2)
            plt.xlabel(r'$ '+vert+' / H $'); 
            plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $')
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
        cs = plt.contourf(X, Y, snap[0,k,:,:]); csName = 'xyCon_' 
        cbar = plt.colorbar()
        plt.xlabel(r'$ x / H $'); 
        plt.ylabel(r'$'+span+' / H $'); 
        plt.tight_layout()
        #fig.show()
        #plt.savefig(figdir+'xy_'+str(k)+'.png', dpi=100)
	mySaveFig(csName, 0)

if snap_plot_yz:
    snap = np.load(datdir+'snap.npy')
    Y, Z = np.meshgrid(y, z)
    for i in range(2,3):
        scale = 3.0;
        fig = plt.figure(figsize=(scale*Ly,scale*Lz))
        cs = plt.contourf(Y, Z, snap[0,:,:,i]);  csName = 'yzCon_'
        #cs = plt.pcolor(Y, Z, snap[0,:,:,i]);  csName = 'yzCol_'
        cbar = plt.colorbar()
        plt.xlabel(r'$'+span+' / H $'); 
        plt.ylabel(r'$'+vert+' / H $'); 
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


if energy_bal:
    dz = Lz / nz
    dx = Lx / nx
    def ddz_w(f):
        # f is on w-nodes, dfdz on uvp-nodes
        dfdz = np.zeros((nx,ny,nz))
        for jz in range(0,nz-1):
            dfdz[jz,:,:] = (f[jz+1,:,:] - f[jz,:,:]) / dz
        return dfdz

    def ddx(f):
        dfdx = np.zeros((nx,ny,nz))
        kx = np.arange(0, nx/2+1)*(1j)
        for j in range(0, ny):
            for k in range(0, nz):
                temp = np.fft.rfft( f[k,j,:] )
                temp = temp * kx
                dfdx[k,j,:] = np.fft.irfft( temp )
        return dfdx

    def ddy(f):
        dfdy = np.zeros((nx,ny,nz))
        ky = np.arange(0, ny/2+1)*(1j)
        for i in range(0, nx):
            for k in range(0, nz):
                temp = np.fft.rfft( f[k,:,i] )
                temp = temp * ky
                dfdy[k,:,i] = np.fft.irfft( temp )
        return dfdy
 
    def interp2uv(f):
        # f is on w-nodes, fuv on uvp-nodes
        fuv = np.zeros((nx,ny,nz))
        for jz in range(0,nz-1):
            fuv[jz,:,:] = (f[jz+1,:,:] + f[jz,:,:]) / 2.0
        return fuv

    def interp2w(f):
        # f is on uv-nodes, fw on w-nodes
        fw = np.zeros((nx,ny,nz))
        for jz in range(0,nz-1):
            fw[jz+1,:,:] = (f[jz,:,:] + f[jz+1,:,:]) / 2.0
        return fw

    dz = Lz / nz
    z_w  = np.linspace(0, Lz, nz, endpoint=False)  # w-grid (no point at z = 1.0)
    z_uv = np.linspace(dz/2, Lz-dz/2, nz) # uv-grid ( last point at Lz-dz/2 )

    u = np.load(datdir+'u.npy')
    v = np.load(datdir+'v.npy')
    w = np.load(datdir+'w.npy')
    txx = np.load(datdir+'txx.npy')
    tyy = np.load(datdir+'tyy.npy')
    tzz = np.load(datdir+'tzz.npy')
    txz = np.load(datdir+'txz.npy')
    txy = np.load(datdir+'txy.npy')
    tyz = np.load(datdir+'tyz.npy')
    rs11 = np.load(datdir+'rs11.npy'); #rs11 = -rs11
    rs22 = np.load(datdir+'rs22.npy'); #rs22 = -rs22
    rs33 = np.load(datdir+'rs33.npy'); #rs33 = -rs33
    rs12 = np.load(datdir+'rs12.npy'); #rs12 = -rs12
    rs13 = np.load(datdir+'rs13.npy'); #rs13 = -rs13
    rs23 = np.load(datdir+'rs23.npy'); #rs23 = -rs23
    rs21 = rs12
    rs31 = rs13
    rs32 = rs23

    dudx = ddx(u)
    dvdx = ddx(v)
    dwdx = ddx(w)
    dudy = ddy(u)
    dvdy = ddy(v)
    dwdy = ddy(w)
    dudz = ddz_w(u); dudz = interp2w(dudz)
    dvdz = ddz_w(v); dvdz = interp2w(dvdz)
    dwdz = ddz_w(w); dwdz = interp2w(dwdz)

    # advection
    E = 0.5 * (u**2 + v**2 + w**2)
    a1 = u*ddx(E)
    b1 = v*ddy(E)
    c1 = ddz_w(E); c1 = w*interp2w(c1)
    term1 = a1 + b1 + c1
    term1 = -1*np.mean(np.mean(term1,axis=2),axis=1)

    # pressure
    dpdx = np.load(datdir+'dpdx.npy')
    dpdx = dpdx - 1.0
    dpdx = interp2w(dpdx)
    dpdy = np.load(datdir+'dpdy.npy')
    dpdy = interp2w(dpdy)
    dpdz = np.load(datdir+'dpdz.npy')
    a2 = u*dpdx
    b2 = v*dpdy
    c2 = w*dpdz
    term2 = a2 + b2 + c2
    term2 = -np.mean(np.mean(term2,axis=2), axis=1)

    fx = np.load(datdir+'fx.npy')
    fx = interp2w(fx)
    fy = np.load(datdir+'fy.npy')
    fy = interp2w(fy)
    fz = np.load(datdir+'fz.npy')

    a5 = u*fx
    b5 = v*fy
    c5 = w*fz
    term5 = a5 + b5 + c5
    term5 = np.mean(np.mean(term5,axis=2),axis=1)

    a6 = u*rs11 + v*rs21 + w*rs31
    b6 = u*rs12 + v*rs22 + w*rs32
    c6 = u*rs13 + v*rs23 + w*rs33
    a6 = ddx(a6)
    b6 = ddy(b6)
    c6 = ddz_w(c6); c6 = interp2w(c6)
    term6a = a6 + b6 + c6
    a6b = u*txx + v*txy + w*txz
    b6b = u*txy + v*tyy + w*tyz
    c6b = u*txz + v*tyz + w*tzz
    a6b = ddx(a6b)
    b6b = ddy(b6b)
    c6b = ddz_w(c6b); c6b = interp2w(c6b)
    term6b = a6b + b6b + c6b 
    term6 = -1*np.mean(np.mean(term6a+term6b,axis=2),axis=1)

    a7 = (rs11+txx) * dudx
    b7 = (rs22+tyy) * dvdy
    c7 = (rs33+tzz) * dwdz
    d7 = (rs12+txy) * dudy
    e7 = (rs13+txz) * dudz
    f7 = (rs23+tyz) * dvdz
    g7 = (rs21+txy) * dvdx
    h7 = (rs31+txz) * dwdx
    i7 = (rs32+tyz) * dwdy
    term7 = a7 + b7 + c7 + d7 + e7 + f7 + g7 + h7 + i7
    term7 = np.mean(np.mean(term7,axis=2),axis=1)

    my_lw = 3    # line width
    my_ms = 5    # marker size
    my_mew = 2   # marker edge width 
    my_fs = 'none'  # fillstyle

    fig = plt.figure()
    plt.plot(term1, z_w, '+b', ms=my_ms, mew=my_mew, fillstyle=my_fs, label=r'$ -\rho \bar{u} \cdot \nabla E $')
    plt.plot(term2, z_w, '^c', ms=my_ms, mew=my_mew, fillstyle=my_fs, label=r'$ -\nabla \cdot (\bar{u}p) $')
    #plt.plot(term3, z_w, 'om', ms=my_ms, mew=my_mew, fillstyle=my_fs, label=r'$ \rho \nu_T \nabla^2 E $')
    #plt.plot(term4, z_w, 'dg', ms=my_ms, mew=my_mew, fillstyle=my_fs, label=r'$ -\rho \nu_T \nabla\bar{u}:\nabla\bar{u} $')
    plt.plot(term5, z_w, '^y', ms=my_ms, mew=my_mew, fillstyle=my_fs, label=r'$ \bar{u} \cdot F $')
    plt.plot(term6, z_w, 'og', ms=my_ms, mew=my_mew, fillstyle=my_fs, label=r'$ \nabla \cdot (\bar{u} \cdot R)+\nabla \cdot (\bar{u}\bar{\tau}) $')
    plt.plot(term7, z_w, 'sr', ms=my_ms, mew=my_mew, fillstyle=my_fs, label=r'$ -R : \nabla \bar{u} + \bar{\tau} \nabla \bar{u} $')

    term_tot = term1 + term2 + term5 + term6 + term7 #+ term3 + term4
    plt.plot(term_tot, z_w, 'ok', ms=my_ms, mew=my_mew, fillstyle=my_fs, label=r'sum')

    xa = -120;   xb = -xa;
    ya = 0;     yb = 1
    line = np.linspace(xa,xb,1000)
    line2 = np.linspace(ya,yb,1000)
    plt.plot(line, np.ones(len(line))*0.15, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.05, '--', color='k')
    plt.plot(np.zeros(len(line2)), line2, '-k', color='k')   
    plt.title('Kinetic energy balance')
    plt.text(xa+10, yb-.1, r'$E = \bar{u} \cdot \bar{u} / 2$',fontsize=24)
    plt.text(xa+10, yb-.2, r'$R = -\rho \overline{u u}$',fontsize=24)
    plt.ylabel(r'$ y / H $')
    plt.xlim((xa,xb)) 
    plt.ylim((ya,yb)) 
    plt.legend(fontsize=18)
    mySaveFig('energy_bal_', 0)


if test_plot:

    def ddx(f, n):
        dfdx = np.fft.rfft(f)
        kx = np.arange(0, n/2+1)*(1j)
        dfdx = dfdx * kx
        dfdx = np.fft.irfft(dfdx)
        return dfdx

    x  = np.linspace(0, Lx, nx, endpoint=False)
    myk = 3.0
    fin = 1.2 + np.sin(myk*x)
    fp1 = myk*np.cos(myk*x)
    fp2 = ddx(fin, nx)

    fig = plt.figure()
    plt.plot(x, fin, 'ok')
    plt.plot(x, fp1, 'sr')
    plt.plot(x, fp2, 'or')
    mySaveFig('test_', 0)



