# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 13:49:27 2020

@author: iveli
"""

import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import pandas
import csv
import os
#from matplotlib.widgets import LassoSelector

#%%
plt.close()
plt.close()
folder_path = r'C:\Users\iveli\OneDrive\Documents\Ron\031919'
file_paths = []

for file in os.listdir(folder_path):
    if file.endswith('.ASC'):
        file_paths.append(folder_path+'\\'+str(file))
        
temps = [1.52, 0.62, 5.71, 3.96, 8, 10, 12, 14, 18, 3.27, 2.24, 0.96];

[val,phase,ffinvspl,chispl,ff,chi] = FFTcalc(fn = file_paths[0])
plt.plot(val)
plt.figure()
plt.plot(1/np.array(ff),chi)
plt.plot(ffinvspl,chispl)


#%%
def FFTcalc(FFrange = [10, 50],FFupi = [], chiupi = [], fn = '', FFmin = 1.5):
    if fn != '':
        try:
            [FFupi, chiupi] = dHvALoad(fn,FFmin)
            
        except:
            print("Failed at reading filename step")
    FFupi = np.array(FFupi)
    FFupi[FFupi < FFrange[0]] = 0
    FFupi[FFupi > FFrange[1]] = 0
    FFup = []
    chiup = []
    for ii in range(len(FFupi)):
        if FFupi[ii] != 0:
            FFup.append(FFupi[ii])
            chiup.append(chiupi[ii])
            
            
    
    n = int(2**11);
    N = int(2**20)
    FFinv = 1/np.array(FFup)
    FFinvspl = np.linspace(min(FFinv),max(FFinv),n) #spline for fitting
    delta_spl = FFinvspl[1] - FFinvspl[0] #int in Ron's code
    
    #For removing curvature in data:
    degree = 3
    chifit = np.polyfit(FFup,chiup,degree) #yfit in ron's code
    chifit_v = np.polyval(chifit,FFup)
    chi = chiup - chifit_v
    
    #Reindexing stuff:
    [toss, unique_ind] = np.unique(FFinv, return_index = True)
    unique_ind = np.sort(unique_ind).astype(int)
    FF = (np.array(FFup))[unique_ind]
    chi = (np.array(chiup))[unique_ind]
    
    chispl_fn = interpolate.interp1d(1/FF, chi)  #np.interp(FFinvspl,1/FF,chi)
    chispl = chispl_fn(FFinvspl);
    n_datapoints = len(chispl) #Datapoints in ron's code
    hwind = np.hanning(n_datapoints)
    chisplW = chispl*hwind; #ysplW in ron's code
    
    #padding stuff with zeros??:
    nZero = int((N-n)/2)
    fft_0 = FFinvspl[0] - nZero*delta_spl #xfft1 in ron' code
    fft_end = FFinvspl[-1] + nZero*delta_spl #xfftEnd in ron's code
    zeros_0 = np.zeros((1,nZero)).flatten() #zeros1 in ron's code
    zeros_end = np.zeros((1,nZero)).flatten() #zerosEnd in ron's code
    ffinv_fft = np.linspace(fft_0,fft_end, N) #xfft in ron's code
    chi_fft = np.concatenate([zeros_0,chisplW,zeros_end]) #yfft in ron's code
    
    n_datapoints = len(chi_fft)
    FFTcomplex = np.fft.fft(chi_fft)
    FFTphase = np.angle(FFTcomplex)
    FFTphase = FFTphase[0:int(n_datapoints/2 + 1)]
    FFTval = abs(FFTcomplex[0:int(n_datapoints/2+1)])
    
    length = abs(ffinv_fft[-1] - ffinv_fft[1])
    fs = n_datapoints/length #sampling frequency
    f = fs/2 * np.linspace(0,1,n_datapoints/2 + 1); #frequency
    
    return FFTval,FFTphase, FFinvspl,chispl,FFup,chiup ;
   
    
    
#%%
def dHvALoad(fn = '', FFmin = 1.5):
    data = pandas.read_csv(fn, delimiter = '\t')
    Bdot = np.array(data.Bdot) #dH/dt
    Mdot = np.array(data.Mag1) #dM/dt (Mag1 in ron's code)
    FF = np.array(data.Field_fixed)
    chi = Mdot/Bdot #(xi in ron's code)
    FFmax = max(FF)
    FFmax_ind = FF.argmax()
    FF[FF < FFmin] = 0
    FFup = FF[0:FFmax_ind]
    chiup = [];
    FFup_return = []
    for ii in range(len(FFup)):
        if FFup[ii] != 0:
            chiup.append(chi[ii])
            FFup_return.append(FFup[ii])
    
    return [FFup_return,chiup]

#%%
def MassCalc():