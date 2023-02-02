# -*- coding: utf-8 -*-
"""
=======================================================
Program : resonator_fluxfit/fluxfit_resonator.py
=======================================================
Summary:
"""
__author__ =  ["James Farmer", "Sadman Ahmed Shanto"]
__date__ = "02/01/2023"
__email__ = "shanto@usc.edu"

#libraries used
import Labber
import numpy as np
import os
import matplotlib.pyplot as plt
from fitTools.Resonator import Resonator
import logging
from scipy.optimize import curve_fit

# flux fitting functions for KO-1
def multiFluxFunc(fl,w0,q0,fl0,fs):
    return w0*(1+q0*np.sin(np.pi*(fl-fl0)/fs/2)*np.arctanh(np.sin(np.pi*(fl-fl0)/fs/2))/(1-np.sin(np.pi*(fl-fl0)/fs/2)*np.arctanh(np.sin(np.pi*(fl-fl0)/fs/2))))**(-0.5)

def FluxFunc(fl,q0):
    return (1+q0*np.sin(np.pi*fl/2)*np.arctanh(np.sin(np.pi*fl/2))/(1-np.sin(np.pi*fl/2)*np.arctanh(np.sin(np.pi*fl/2))))**(-0.5)

if __name__ == "__main__":
    plt.rcParams.update({'font.size':14})

    fpath = input("Path to .hdf5 file: ").replace('"','')

    path,fname = os.path.split(fpath)
    path += r'\\'
    figpath = 'figures\\'+fname[:-4]+'\\'
    if not os.path.exists(path+'figures\\'):
        os.mkdir(path+'figures\\')
    if not os.path.exists(path+figpath):
        os.mkdir(path+figpath)
    lf = Labber.LogFile(path + fname)

    logFileName = path + f"fluxfit_info_{fname[:-4]}log"
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler(logFileName, 'a'))
    print = logger.info

    nEntries = lf.getNumberOfEntries()
    try:
        current = np.squeeze(np.round(lf.getData(name='Victor - Source current'),decimals=6))*1e3
    except:
        try:
            current = np.squeeze(np.round(lf.getData(name='Vikram - Source current'),decimals=6))*1e3
        except:
            current = np.squeeze(np.round(lf.getData(name='Vladimir - Source current'),decimals=6))*1e3



    fits = {'f':[],'Q':[],'Qint':[],'Qext':[]}
    for n in range(nEntries):
        (xdata,ydata) = lf.getTraceXY(entry = n)
        phase = np.unwrap(np.angle(ydata))
        dphase = np.diff(phase)
        f0_est_ind = np.argmax(np.abs(dphase))
        f0_est = xdata[f0_est_ind]
        normphase = phase - phase[f0_est_ind]
        segment_mask = np.logical_and(normphase > -np.pi/2, normphase < np.pi/2)
        kappa_est = xdata[segment_mask][-1] - xdata[segment_mask][0]
        fstart = (f0_est - 2*kappa_est)
        fstop = (f0_est + 2*kappa_est)
        mask = np.logical_and(xdata >= fstart, xdata <= fstop)
        res = Resonator('r',xdata[mask],ydata[mask])
        res.autofit(electric_delay=0)
        fits['f'].append(res.f0*1e-9)
        fits['Q'].append(res.Q)
        fits['Qint'].append(res.Qi)
        fits['Qext'].append(res.Qc)

        if res.fit_found:
            print(20*"=")
            res.show(savefile = path+figpath+fname[:-4]+'resonance_{}-mA.png'.format(current[n]))
            print('\nFit at {} mA'.format(current[n]))
            print(res)
            print(20*"="+"\n")
            
    # fit the flux curves
    pars,covs = curve_fit(multiFluxFunc,current[cmask],freqs[cmask],p0=[5.72,0.01,-13,20])

    fig = plt.figure(figsize=[9,6],constrained_layout=True)
    plt.plot(current,fits['f'],'r.')
    plt.title('Frequency vs current')
    plt.xlabel('Source current [mA]')
    #plt.xlabel('LO attenuation [dB]')
    plt.ylabel('Frequency from fit [GHz]')
    plt.savefig(path+figpath+fname[:-4]+'_f0-vs-current.png')
    plt.show()

    fig = plt.figure(figsize=[9,6],constrained_layout=True)
    plt.scatter(current,fits['Q'],s=20,c='r',label='Total Q')
    plt.scatter(current,fits['Qint'],s=20,c='b',label='Internal Q')
    plt.scatter(current,fits['Qext'],s=20,c='g',label='external Q')
    plt.title('Q vs current')
    plt.xlabel('Source current [mA]')
    #plt.xlabel('LO attenuation [dB]')
    plt.ylabel('Quality factor')
    plt.legend()
    plt.savefig(path+figpath+fname[:-4]+'_Q-vs-Current.png')
    plt.show()


    print(f"\n\nLog File: {logFileName}")
    print(f"\nPlots Directory: {path + figpath}")







