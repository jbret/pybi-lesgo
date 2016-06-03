"""
Python plotting for lesgo binary data.
Author: Joel Bretheim
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
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

Lx = 2*np.pi;
Ly = 2*np.pi;
Lz = 1.0;
#kx_vec=[0,4,8]
kx_vec=[0]
#kx_vec=[0,1,2]

vel_avg_plot = 1;
rs_plot      = 1;
sp_plot      = 1;
sp1d_plot    = 1;
vel2_rs_plot = 1;
tau_plot     = 1;
spanSpec_plot= 1;
snap_plot    = 0;  thisSnap = 5000;  # on uv-grid
snap_plot_yz = 0;
snap_plot_xy = 0;
fourier      = 1;

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
    #plt.semilogx(z, uMean, 'o')
    plt.semilogx(z/(1./180), uMean, 'o')
    #plt.semilogx(z, 1/0.4*np.log(z/.0001), '-k', label=r'$1/\kappa \ \mathrm{log}(z/z_{0})$')
    plt.semilogx(z/(1./180), 1/0.41*np.log(z/(1./180))+5.0, '-k', label=r'$1/\kappa \ \mathrm{log}(z^{+})+B$')
    #plt.xlabel('$ z / H $', fontsize=18);
    plt.xlabel('$z^{+} $', fontsize=18);
    plt.ylabel('$[\ u / u_*]$', fontsize=18); 
    #plt.xlim([.02 ,1.1])
    plt.xlim([1, 1000])
    #plt.text(.4,3,r'$ \kappa = 0.4,\ z_{0} = 10^{-4} $', fontsize=14)
    plt.text(110,2.5,r'$ \kappa = 0.41,\ B = 5.0 $', fontsize=14)
    plt.legend(loc='lower right', fontsize=14)
    plt.tight_layout()
    mySaveFig('mvp_', 0)

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
    plt.plot(-1*(rs13Mean + txzMean), z, '-s', color='k', markeredgewidth=1, markerfacecolor="None", label = r'$ \mathrm{sum} $')
    line = np.linspace(0,1,1000)
    plt.plot(line, 1+-1.0*line, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.15, '--', color='k')
    plt.plot(line, np.ones(len(line))*0.05, '--', color='k')
    plt.xlabel(r'$ \mathrm{Stress} $', fontsize=18); plt.ylabel(r'$ z / H $', fontsize=18)
    plt.legend()
    plt.tight_layout()
    mySaveFig('tau_', 0)


if spanSpec_plot:
    sp11_1d = np.load(datdir+'sp11_1d.npy')
    sp22_1d = np.load(datdir+'sp22_1d.npy')
    sp33_1d = np.load(datdir+'sp33_1d.npy')

    spec11 = np.mean(sp11_1d[:,:,:], axis=2)
    spec22 = np.mean(sp22_1d[:,:,:], axis=2)
    spec33 = np.mean(sp33_1d[:,:,:], axis=2)

    #spec11 = np.mean(sp11_1d[:,:,:], axis=2)
    #spec11 = np.sum(spec11[:,0:], axis=1)
    #spec22 = np.mean(sp22_1d[:,:,:], axis=2)
    #spec22 = np.sum(spec22[:,0:], axis=1)
    #spec33 = np.mean(sp33_1d[:,:,:], axis=2)
    #spec33 = np.sum(spec33[:,0:], axis=1)

    heights = [ 5, 10, 30 ];  numH = np.size(heights)
    ky = np.arange(0,ny/2)
    fig = plt.figure(figsize=(10,8))
    for j in range(0,numH):
        plt.subplot(1,numH,j+1)
        k = heights[j]
        plt.loglog(ky[1:], spec11[k,1:ny/2],'o', label=r'$E_{uu}$')
        plt.loglog(ky[1:], spec22[k,1:ny/2],'o', label=r'$E_{vv}$')
        plt.loglog(ky[1:], spec33[k,1:ny/2],'o', label=r'$E_{ww}$')
        plt.title(r'$ z^{+} = $'+str(heights[j]));
        plt.xlabel(r'$ k_y $'); plt.tight_layout()
        plt.legend(loc='lower left')
        #plt.ylim([])

    #plt.plot(z[1:], comp11[1:], 's')
    mySaveFig('ky_spec', 0)


if sp_plot:
    sp11 = np.load(datdir+'sp11.npy')
    sp22 = np.load(datdir+'sp22.npy')
    sp33 = np.load(datdir+'sp33.npy')
    ky = np.arange(0,ny/2)
    lamY = Ly / ky

    fig = plt.figure()
    plt.loglog(ky, sp11[5,0:ny/2,0],'o', color='g',label = r'$ E_{uu}, z^{+}=5 $')
    plt.loglog(ky, sp22[5,0:ny/2,0],'o', color='b',label = r'$ E_{vv}, z^{+}=5 $')
    plt.loglog(ky, sp33[5,0:ny/2,0],'o', color='r',label = r'$ E_{ww}, z^{+}=5 $')
    plt.loglog(ky, sp11[19,0:ny/2,0],'s', color='g',label = r'$ E_{uu}, z^{+}=19 $')
    plt.loglog(ky, sp22[19,0:ny/2,0],'s', color='b',label = r'$ E_{vv}, z^{+}=19 $')
    plt.loglog(ky, sp33[19,0:ny/2,0],'s', color='r',label = r'$ E_{ww}, z^{+}=19 $')
    plt.loglog(ky, sp11[100,0:ny/2,0],'^', color='g',label = r'$ E_{uu}, z^{+}=100 $')
    plt.loglog(ky, sp22[100,0:ny/2,0],'^', color='b',label = r'$ E_{vv}, z^{+}=100 $')
    plt.loglog(ky, sp33[100,0:ny/2,0],'^', color='r',label = r'$ E_{ww}, z^{+}=100 $')
    plt.xlabel(r'$ k_{y} $'); plt.ylabel(r'$ E_{ii} $')
    plt.ylim([10**(-14), 10**1])
    plt.legend(loc='lower left',ncol=3)
    plt.tight_layout()
    mySaveFig('spanSpec_', 0)

    #levels=[0.0,0.05,0.1,0.15,0.2,0.25,0.3]
    sp11sum1 = np.sum(sp11[:,:,1:],axis=2)  # sum over kx
    #sp11sum13 = np.sum(sp11sum1,axis=0) # sum over z
    #kyE11 = ky * sp11sum1[:,0:ny/2]/ sp11sum13[0:ny/2]
    kyE11 = ky * sp11sum1[:,0:ny/2]

    sp22sum1 = np.sum(sp22[:,:,1:],axis=2)  # sum over kx
    kyE22 = ky * sp22sum1[:,0:ny/2]

    sp33sum1 = np.sum(sp33[:,:,1:],axis=2)  # sum over kx
    kyE33 = ky * sp33sum1[:,0:ny/2]

    #kyE33 = np.zeros([nz,ny/2])
    #for k in range(0,nz):
    #    for j in range (0,ny/2):
    #        kyE33[k,j] = ky[j]*sp22sum1[k,j]/sp22sum13[j]

    LAMY, Z = np.meshgrid(lamY[1:], z[1:])

    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(1,3,1)
    cs = plt.contourf(LAMY, Z, kyE11[1:,1:]);
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ \lambda_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$k_y E_{uu}$')
    plt.tight_layout()
     
    ax = fig.add_subplot(1,3,2)
    cs = plt.contourf(LAMY, Z, kyE22[1:,1:]);  
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ \lambda_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$k_y E_{vv}$')
    plt.tight_layout()
    
    ax = fig.add_subplot(1,3,3)
    cs = plt.contourf(LAMY, Z, kyE33[1:,1:]);  
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ \lambda_y / \delta $', fontsize=18);
    plt.ylabel(r'$ z / \delta $', fontsize=18); 
    plt.title(r'$k_y E_{ww}$')
    plt.tight_layout()
    mySaveFig('kyEii_', 0)

    lamYp = lamY/(1./180)
    zp = z/(1./180)
    LAMYP, ZP = np.meshgrid(lamYp[1:], zp[1:])
    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(1,3,1)
    cs = plt.contourf(LAMYP, ZP, kyE11[1:,1:]);
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ \lambda_y^{+} $', fontsize=18);
    plt.ylabel(r'$ z^{+} $', fontsize=18); 
    plt.title(r'$k_y E_{uu}$')
    plt.tight_layout()
     
    ax = fig.add_subplot(1,3,2)
    cs = plt.contourf(LAMYP, ZP, kyE22[1:,1:]);  
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ \lambda_y^{+} $', fontsize=18);
    plt.ylabel(r'$ z^{+} $', fontsize=18); 
    plt.title(r'$k_y E_{vv}$')
    plt.tight_layout()
    
    ax = fig.add_subplot(1,3,3)
    cs = plt.contourf(LAMYP, ZP, kyE33[1:,1:]);  
    ax.set_xscale('log');  ax.set_yscale('log')
    cbar = plt.colorbar()
    plt.xlabel(r'$ \lambda_y^{+} $', fontsize=18);
    plt.ylabel(r'$ z^{+} $', fontsize=18); 
    plt.title(r'$k_y E_{ww}$')
    plt.tight_layout()
    mySaveFig('kyEiiplus_', 0)
    
if rs_plot:
    rs11Mean = np.load(datdir+'rs11Mean.npy')
    rs22Mean = np.load(datdir+'rs22Mean.npy')
    rs33Mean = np.load(datdir+'rs33Mean.npy')
    rs13Mean = np.load(datdir+'rs13Mean.npy')
    rs23Mean = np.load(datdir+'rs23Mean.npy')
    rs12Mean = np.load(datdir+'rs12Mean.npy')
    fig = plt.figure()
    plt.subplot(2,1,1)
    L1=plt.plot(z, rs11Mean, 'o', color='b', label = r'$[ u^{\prime} u^{\prime}]$')
    L2=plt.plot(z, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    L3=plt.plot(z, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    L4=plt.plot(z, rs13Mean, 'o', color='c', label = r'$[ u^{\prime} w^{\prime}]$')
    L5=plt.plot(z, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    L6=plt.plot(z, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    plt.legend(loc='upper right',ncol=2)
    plt.xlabel(r'$ z / H $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.tight_layout()
    plt.subplot(2,1,2)
    Re = 180; zp = z/(1./Re)
    L1=plt.plot(zp, rs11Mean, 'o', color='b', label = r'$[ u^{\prime} u^{\prime}]$')
    L2=plt.plot(zp, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    L3=plt.plot(zp, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    L4=plt.plot(zp, rs13Mean, 'o', color='c', label = r'$[ u^{\prime} w^{\prime}]$')
    L5=plt.plot(zp, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    L6=plt.plot(zp, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    plt.legend(loc='upper right',ncol=2)
    plt.xlabel(r'$ z^{+} $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.tight_layout()
    mySaveFig('rs_', 0)

    fig = plt.figure()
    plt.subplot(2,1,1)
    L2=plt.plot(z, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    L3=plt.plot(z, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    L4=plt.plot(z, -1.0*rs13Mean, 'o', color='c', label = r'$[ -u^{\prime} w^{\prime}]$')
    L5=plt.plot(z, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    L6=plt.plot(z, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    plt.legend(loc='upper center',ncol=5)
    plt.xlabel(r'$ z / H $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.tight_layout()
    plt.subplot(2,1,2)
    Re = 180; zp = z/(1./Re)
    L2=plt.plot(zp, rs22Mean, 'o', color='g', label = r'$[ v^{\prime} v^{\prime}]$')
    L3=plt.plot(zp, rs33Mean, 'o', color='r', label = r'$[ w^{\prime} w^{\prime}]$')
    L4=plt.plot(zp, -1.0*rs13Mean, 'o', color='c', label = r'$[ -u^{\prime} w^{\prime}]$')
    L5=plt.plot(zp, rs23Mean, 'o', color='m', label = r'$[ v^{\prime} w^{\prime}]$')
    L6=plt.plot(zp, rs12Mean, 'o', color='k', label = r'$[ u^{\prime} v^{\prime}]$')
    plt.legend(loc='upper center',ncol=5)
    plt.xlabel(r'$ z^{+} $', fontsize=18); plt.ylabel(r'$ [ u_{i}^{\prime} u_{i}^{\prime}]/u_{*}^{2} $', fontsize=18)
    plt.tight_layout()
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


            #sp11_1d = np.load(datdir+'sp11_1d.npy')
            #fig = plt.figure()
            #plt.plot(z,rs11Mean,'o')
            #plt.plot(z,uuMean,'^')
            #comp11 = np.mean(sp11_1d[:,:,:], axis=2)
            #comp11 = np.sum(comp11[:,0:], axis=1)
            #plt.plot(z[1:], comp11[1:], 's')
            #mySaveFig('comp', 0)


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
