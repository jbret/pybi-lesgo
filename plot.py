"""
Python plotting for lesgo binary data.
Author: Joel Bretheim
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

Lx = np.pi;
Ly = np.pi;
Lz = 1.0;
#kx_vec=[0,4,8]
kx_vec=[0]
#kx_vec=[0,1,2]

vel_avg_plot = 1;
rs_plot      = 1;
tau_plot     = 1;
snap_plot    = 0;  thisSnap = 5000;  # on uv-grid
snap_plot_yz = 0;
snap_plot_xy = 0;
fourier      = 1;

z = np.linspace(0, Lz, nz, endpoint=True)
y = np.linspace(0, Ly, ny, endpoint=False)
x = np.linspace(0, Lx, nx, endpoint=False)

datdir = 'pyData/'
figdir1 = 'pyFigs/'
figdir2 = figdir1 + 'otherFormats/'
figdir3 = figdir1 + 'highDPI/'
system('mkdir ' + figdir1)
system('mkdir ' + figdir2)
system('mkdir ' + figdir3)

plt.close("all")
# begin plotting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if vel_avg_plot:
    uMean = np.load(datdir+'uMean.npy')
    fig = plt.figure()
    plt.semilogx(z, uMean, 'o')
    plt.semilogx(z, 1/0.4*np.log(z/.0001), '-k', label=r'$1/\kappa \ \mathrm{log}(z/z_{0})$')
    plt.xlabel('$z / H$', fontsize=18); plt.ylabel('$[\ u / u_*]$', fontsize=18); 
    plt.xlim([.02 ,1.1])
    plt.text(.4,3,r'$ \kappa = 0.4,\ z_{0} = 10^{-4} $', fontsize=14)
    plt.legend(loc='lower right', fontsize=14)
    plt.tight_layout()
    plt.savefig(figdir1+'mvp_' + runName + '.png')
    plt.savefig(figdir2+'mvp_' + runName + '.pdf')
    plt.savefig(figdir2+'mvp_' + runName + '.eps')
    plt.savefig(figdir2+'mvp_' + runName + '.jpg')
    plt.savefig(figdir3+'dpi600-'+'mvp_' + runName + '.png',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'mvp_' + runName + '.pdf',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'mvp_' + runName + '.eps',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'mvp_' + runName + '.jpg',dpi=600)
    #fig.show()

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
    ax.add_artist(circle1)
    ax.add_artist(circle2)
    ax.add_artist(circle3)
    ax.add_artist(circle4)
    ax.add_artist(circle5)
    ax.add_artist(circle6)
    plt.tight_layout()
    plt.savefig(figdir1+'uXmean_' + runName + '.png')
    plt.savefig(figdir2+'uXmean_' + runName + '.pdf')
    plt.savefig(figdir2+'uXmean_' + runName + '.eps')
    plt.savefig(figdir2+'uXmean_' + runName + '.jpg')
    plt.savefig(figdir3+'dpi600-'+'uXmean_' + runName + '.png',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'uXmean_' + runName + '.pdf',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'uXmean_' + runName + '.eps',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'uXmean_' + runName + '.jpg',dpi=600)
    #fig.show()

if tau_plot:
    rs13Mean = np.load(datdir+'rs13Mean.npy')
    txzMean = np.load(datdir+'txzMean.npy')
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
    plt.savefig(figdir1+'tau_' + runName + '.png')
    plt.savefig(figdir2+'tau_' + runName + '.pdf')
    plt.savefig(figdir2+'tau_' + runName + '.eps')
    plt.savefig(figdir2+'tau_' + runName + '.jpg')
    plt.savefig(figdir3+'dpi600-'+'tau_' + runName + '.png',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'tau_' + runName + '.pdf',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'tau_' + runName + '.eps',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'tau_' + runName + '.jpg',dpi=600)
    #fig.show()
    
if rs_plot:
    rs11Mean = np.load(datdir+'rs11Mean.npy')
    rs22Mean = np.load(datdir+'rs22Mean.npy')
    rs33Mean = np.load(datdir+'rs33Mean.npy')
    rs13Mean = np.load(datdir+'rs13Mean.npy')
    rs23Mean = np.load(datdir+'rs23Mean.npy')
    rs12Mean = np.load(datdir+'rs12Mean.npy')
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
    plt.savefig(figdir1+'rs_' + runName + '.png')
    plt.savefig(figdir2+'rs_' + runName + '.pdf')
    plt.savefig(figdir2+'rs_' + runName + '.eps')
    plt.savefig(figdir2+'rs_' + runName + '.jpg')
    plt.savefig(figdir3+'dpi600-'+'rs_' + runName + '.png',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'rs_' + runName + '.pdf',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'rs_' + runName + '.eps',dpi=600)
    plt.savefig(figdir3+'dpi600-'+'rs_' + runName + '.jpg',dpi=600)
    #fig.show()

if snap_plot_xy:
    snap = np.load(datadir+'snap.npy')
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
    snap = np.load(datadir+'snap.npy')
    Y, Z = np.meshgrid(y, z)
    for i in range(2,3):
        scale = 3.0;
        fig = plt.figure(figsize=(scale*Ly,scale*Lz))
        cs = plt.contourf(Y, Z, snap[0,:,:,i]);  csName = 'yzCon'
        #cs = plt.pcolor(Y, Z, snap[0,:,:,i]);  csName = 'yzCol'
        cbar = plt.colorbar()
        plt.xlabel(r'$ y / H $', fontsize=18); plt.ylabel(r'$ z / H $', fontsize=18); 
        plt.tight_layout()
        plt.savefig(figdir1 + csName + runName + '.png')
        plt.savefig(figdir2 + csName + runName + '.pdf')
        plt.savefig(figdir2 + csName + runName + '.eps')
        plt.savefig(figdir2 + csName + runName + '.jpg')
        plt.savefig(figdir3+'dpi600-'+ csName + runName + '.png',dpi=600)
        plt.savefig(figdir3+'dpi600-'+ csName + runName + '.pdf',dpi=600)
        plt.savefig(figdir3+'dpi600-'+ csName + runName + '.eps',dpi=600)
        plt.savefig(figdir3+'dpi600-'+ csName + runName + '.jpg',dpi=600)
        #fig.show()
        
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

        

