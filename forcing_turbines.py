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

single_turbine_plot = 0;  turbNum = 1
all_turbine_power   = 1;  total_turb_num = 24;
# turbine numbering progresses along rows first, then columns

# col 0: current time
# col 1: instantaneous disk-averaged velocity
# col 2: current time- and disk-averaged velocity
# col 3: instantaneous total force for this turbine
# col 4: instantaneous power for this turbine

if all_turbine_power:
    plt.close("all")
    t   = np.loadtxt('./turbine/turbine_1_forcing.dat', usecols=(0,) )
    numT = np.size(t)
    pi_all = np.zeros((numT, total_turb_num))
    for i in range(1, total_turb_num + 1 ):
        filename = './turbine/turbine_'+str(i)+'_forcing.dat'
        pi_all[:,i-1] = np.loadtxt(filename, usecols=(4,) )

    farm_avg = np.mean( pi_all[:,:], axis=1)
    t1 = int(np.floor(1./8*numT))
    t2 = numT - t1

    farm_t_avg = np.mean(farm_avg[ t1: ] )

    fig = plt.figure(figsize=(10,5))
    plt.plot(t, farm_avg)
    plt.plot(t[t1:], np.ones(t2)*farm_t_avg,'-r')
    plt.title('farm time avg power = '+str(farm_t_avg))
    plt.ylim([0, 8])
    plt.savefig(runName+'_farm_avg.png')

    #fig = plt.figure(figsize=(10,5))
    #for i in range(1, total_turb_num + 1 ):
    #    plt.plot(t, pi_all[:,i-1])
    #plt.savefig('pi_all.png')

if single_turbine_plot:
    plt.close("all")
    filename1 = './turbine/turbine_'+str(turbNum)+'_forcing.dat'
    if isfile(filename1):
        t   = np.loadtxt(filename1, usecols=(0,) )
        vi1 = np.loadtxt(filename1, usecols=(1,) ); vi1_avg = np.mean(vi1)
        vi2 = np.loadtxt(filename1, usecols=(2,) ); vi2_avg = np.mean(vi2)
        fi  = np.loadtxt(filename1, usecols=(3,) ); fi_avg  = np.mean(fi)
        pi  = np.loadtxt(filename1, usecols=(4,) ); pi_avg  = np.mean(pi)

        fig = plt.figure(figsize=(10,15))
        plt.subplot(411)
        plt.plot(t, vi1, 'k')
        plt.plot(t, np.ones(np.size(t))*vi1_avg,'-r')
        plt.subplot(412)
        plt.plot(t, vi2, 'k')
        plt.plot(t, np.ones(np.size(t))*vi2_avg,'-r')
        plt.subplot(413)
        plt.plot(t, fi, 'k')
        plt.plot(t, np.ones(np.size(t))*fi_avg,'-r')
        plt.subplot(414)
        plt.plot(t, pi, 'k')
        plt.plot(t, np.ones(np.size(t))*pi_avg,'-r')

        plt.xlabel('time'); 
        plt.tight_layout()
        saveString = 'single_turbine_'+str(turbNum) + '_' + runName 
        plt.savefig(saveString + '.png')
    else:
        print ">>>> File "+filename1+" is not present!"

