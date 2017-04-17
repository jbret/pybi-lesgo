"""
Python plotting for lesgo binary data
Author: Joel Bretheim
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
#matplotlib.rc('text', usetex = True) # disable on MARCC (anaconda-python/2.7.10)
import matplotlib.pyplot as plt
from cycler import cycler
from os import getcwd
from os.path import isfile

myDir = getcwd(); dirParts = myDir.split("/")
runName = dirParts[len(dirParts)-1]; print "This run's name: ", runName

filename1  = './output/tau_wall.dat'
filename1a = './output/tau_wall_bot.dat'
filename1b = './output/tau_wall_top.dat'
filename2 = './output/ke_kx.dat'
filename3 = './output/ke_ky.dat'
filename4 = './output/check_ke.dat'

def get_wall_stress(myFileName):
    lines = np.loadtxt(myFileName, skiprows=1, usecols=(0,1,7))
    ts = lines[:,0]        # timesteps
    t = lines[:,1]         # nondimensional time
    tau_wall = lines[:,2]  # nondimensional wallstress
    return ts, t, tau_wall

plt.close("all")
if isfile(filename1):
    ts, t, tau_wall = get_wall_stress(filename1)
if isfile(filename1a):
    ts, t, tau_wall = get_wall_stress(filename1a)

if isfile(filename1) or isfile(filename1a):
    n = np.size(t)
    meanStress = np.mean(tau_wall[np.floor(n/2):n-1]) 
    baseline = np.ones(n) * meanStress
    fig = plt.figure(figsize=(10,7))
    ymin = np.min(tau_wall) - 0.1*np.min(tau_wall)
    ymax = np.max(tau_wall) + 0.1*np.max(tau_wall)
    plt.subplot(2,1,1)
    plt.plot(t, baseline, 'k')
    plt.plot(t, tau_wall, 'or')
    plt.xlabel('nondimensional time'); plt.ylabel('wall shear stress')
    plt.xlim((0, t[-1]))
    plt.ylim((ymin, ymax))
    plt.text(0.0, ymax-0.05*ymax, 'mean stress = '+str(meanStress) ) 
    plt.subplot(2,1,2) 
    plt.plot(ts, baseline, 'k')
    plt.plot(ts, tau_wall, 'or')
    plt.xlabel('timestep'); plt.ylabel('wall shear stress')
    plt.ylim((ymin, ymax))
    plt.tight_layout()
    plt.savefig('hist_tau_wall_bot_' + runName + '.png')
    #fig.show()

if isfile(filename1b):
    ts, t, tau_wall = get_wall_stress(filename1b)
    n = np.size(t)
    meanStress = np.mean(tau_wall[np.floor(n/2):n-1]) 
    baseline = np.ones(n) * meanStress
    fig = plt.figure(figsize=(10,7))
    ymin = np.min(tau_wall) - 0.1*np.min(tau_wall)
    ymax = np.max(tau_wall) + 0.1*np.max(tau_wall)
    plt.subplot(2,1,1)
    plt.plot(t, baseline, 'k')
    plt.plot(t, tau_wall, 'ob')
    plt.xlabel('nondimensional time'); plt.ylabel('wall shear stress')
    plt.text(0.0, ymax-0.05*ymax, 'mean stress = '+str(meanStress) ) 
    plt.xlim((0, t[-1]))
    plt.ylim((ymin, ymax))
    plt.subplot(2,1,2) 
    plt.plot(ts, baseline, 'k')
    plt.plot(ts, tau_wall, 'ob')
    plt.xlabel('timestep'); plt.ylabel('wall shear stress')
    plt.ylim((ymin, ymax))
    plt.tight_layout()
    plt.savefig('hist_tau_wall_top_' + runName + '.png')

if isfile(filename2):
    lines = np.loadtxt(filename2, skiprows=1)
    t = lines[:,0]
    num_kx = np.size(lines,1)
    max_kx = 70

    if num_kx < max_kx:
	kx_plot = num_kx
    else:
  	kx_plot = max_kx

    fig = plt.figure(figsize=(12,6))
    for i in range(1, kx_plot-1):
        plt.semilogy(t, lines[:,i], label=r'$k_x =$'+str(i-1))
        plt.rc('lines', linewidth=2)
        plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
                                   cycler('marker', ['o', 's','o','s'])))
                                   #cycler('linestyle', ['-', '--', ':', '-.'])))

    h = np.size(t)/6
    h2 = np.size(t)-h
    kxavg = np.zeros((h2,kx_plot))
    for i in range(1,kx_plot):
	kxavg[:,i-1] = np.sum(lines[h:,i]) / h2
        plt.semilogy(t[h:], kxavg[:,i-1], '-')

    tot_kx = np.sum(kxavg[0,:])
    tot = np.zeros((h2,1))
    tot[:,0] = tot_kx
    plt.semilogy(t[h:], tot[:,0], '-')

    kxnorm = kxavg / tot_kx
    
    plt.xlabel('timestep'); plt.ylabel('kx energy')
    #plt.legend(loc='lower left',fancybox=True, shadow=True,ncol=10)
    #plt.legend(loc='lower center', fancybox=True, shadow=True, ncol=num_kx-1)
    #lgd = plt.legend(bbox_to_anchor=(0.,1.02,1.,.75), loc=2, ncol=6, mode="expand", borderaxespad=0.)
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=4, borderaxespad=0.)
    ymin, ymax = plt.ylim()
    plt.ylim((10**1, ymax))
    #plt.tight_layout()
    plt.savefig('hist_kx_energy' + runName + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    #fig.show()

    fig = plt.figure(figsize=(12,6))
    #kxs = np.array(range(0,32))
    #plt.plot(kxs,kxnorm[0,0:32],'-og')
    #f = 0.14*kxs**(-1.0)
    #plt.plot(kxs,f,'-ob')
    #plt.savefig('kxs_' + runName + '.png')

else:
    print ">>>> File "+filename2+" is not present!"

if isfile(filename3):
    lines = np.loadtxt(filename3, skiprows=1)
    t = lines[:,0]
    num_ky = np.size(lines,1)
    max_ky = 70

    if num_ky < max_ky:
	ky_plot = num_ky
    else:
  	ky_plot = max_ky

    fig = plt.figure(figsize=(12,6))
    for i in range(1, ky_plot-1):
        plt.semilogy(t, lines[:,i], label=r'$k_y =$'+str(i-1))
        plt.rc('lines', linewidth=2)
        plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
                                   cycler('marker', ['o', 's','o','s'])))
    
    plt.xlabel('timestep'); plt.ylabel('ky energy')
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=4, borderaxespad=0.)
    ymin, ymax = plt.ylim()
    plt.ylim((10**1, ymax))
    plt.savefig('hist_ky_energy' + runName + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
else:
    print ">>>> File "+filename3+" is not present!"

if isfile(filename4):
    lines = np.loadtxt(filename4, skiprows=1, usecols=(0,1))
    t = lines[:,0]        # time
    ke = lines[:,1]
    fig = plt.figure(figsize=(10,7))
    plt.plot(t, ke, 'og')
    plt.ylabel('KE');  plt.xlabel('time');
    plt.xlim((0, t[-1]))
    plt.tight_layout()
    plt.savefig('hist_check_ke_+' + runName + '.png')


