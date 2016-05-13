"""
Python plotting for lesgo binary data
Author: Joel Bretheim, jbretheim@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from os import getcwd
from os.path import isfile

myDir = getcwd(); dirParts = myDir.split("/")
runName = dirParts[len(dirParts)-1]; print "This run's name: ", runName

filename1 = './output/tau_wall.dat'
filename2 = './output/ke_kx.dat'

plt.close("all")
if isfile(filename1):
    lines = np.loadtxt(filename1, skiprows=1, usecols=(0,6,7))
    t = lines[:,0]
    one = lines[:,1]
    tau_wall = lines[:,2]

    fig = plt.figure(figsize=(12,3))
    plt.plot(t, one, 'k')
    plt.plot(t, tau_wall, 'or')
    plt.xlabel('timestep'); plt.ylabel('wall shear stress')
    plt.tight_layout()
    plt.savefig('hist_tau_wall' + runName + '.png')
    fig.show()
else:
    print ">>>> File "+filename1+" is not present!"

if isfile(filename2):
    lines = np.loadtxt(filename2, skiprows=1)
    t = lines[:,0]
    num_kx = np.size(lines,1)

    fig = plt.figure(figsize=(12,3))
    for i in range(1, num_kx-1):
        plt.semilogy(t, lines[:,i], label=r'$k_x =$'+str(i-1))
        plt.rc('lines', linewidth=2)
        plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
                                   cycler('marker', ['o', 's','o','s'])))
                                   #cycler('linestyle', ['-', '--', ':', '-.'])))

    plt.xlabel('timestep'); plt.ylabel('kx energy')
    plt.legend(loc='lower center', fancybox=True, shadow=True, ncol=num_kx-1)
    plt.tight_layout()
    plt.savefig('hist_kx_energy' + runName + '.png')
    fig.show()
else:
    print ">>>> File "+filename2+" is not present!"

