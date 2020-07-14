
import numpy as np
import scipy as sp
from scipy import interpolate,signal, optimize
import matplotlib.pyplot as plt
import pandas
import csv
import os
from enum import Enum
#%%
class ClassdHvALoad():
    FFmin = 1.5
    chiup = []
    FFup = []

class ClassFFTCalc():
    FFrange = [10,50]
    ffInv = []
    chi = []
    ffInvSpline = []
    chiSpline = []
    frequencies= []
    FFT = []
    phase = []
    
class ClassMassCalc():
    m = []
    A = []
    AoT = []
    AoTfitV = []
    tempVec = []
    maxFrequency = []
    

class dataInterp(): 
    inputFrequencyRange = []
    fn = ''
    inputTemperature = []
    dHvALoadObj = ClassdHvALoad()
    FFTCalcObj = ClassFFTCalc()
    MassCalcObj = ClassMassCalc()
    
    def func_dHvALoad(self):
        FFmin = self.dHvALoadObj.FFmin
        data = pandas.read_csv(self.fn, delimiter = '\t')
        Bdot = np.array(data.Bdot) #dH/dt
        Mdot = np.array(data.Mag1) #dM/dt (Mag1 in ron's code)
        FF = np.array(data.Field_fixed)
        chi = Mdot/Bdot #(xi in ron's code)
        #FFmax = max(FF)
        FFmax_ind = FF.argmax()
        FF[FF < FFmin] = 0
        FFup = FF[0:FFmax_ind]
        chiup = [];
        FFup_return = []
        for ii in range(len(FFup)):
            if FFup[ii] != 0:
                chiup.append(chi[ii])
                FFup_return.append(FFup[ii])
        self.dHvALoadObj.chiup = np.array(chiup)
        self.dHvALoadObj.FFup = np.array(FFup)
        #print('Finished loading dHvA data\n')
        
    def func_FFTCalc(self):
        if self.dHvALoadObj.FFup != [] & self.dHvALoadObj.chiup != []:
            FFupi = self.dHvALoadObj.FFup
            chiupi = self.dHvALoadObj.chiup
        else:
            print('Loading dHvA data first\n')
            self.func_dHvALoad(self)
            FFupi = self.dHvALoadObj.FFup
            chiupi = self.dHvALoadObj.chiup
            print('Done loading dHvA data')
        FFupi[FFupi < self.FFTCalcObj.FFrange[0]] = 0
        FFupi[FFupi > self.FFTCalcObj.FFrange[1]] = 0
        nonzeroIndices = np.where(FFupi != 0)[1]
        FFup = FFupi[nonzeroIndices]
        chiup = chiupi[nonzeroIndices]
        
        n = int(2**11);
        N = int(2**20)
        FFinv = 1/FFup
        FFinvspl = np.linspace(min(FFinv),max(FFinv),n) #spline for fitting
        delta_spl = FFinvspl[1] - FFinvspl[0] #int in Ron's code
        
        degree = 3
        chifit = np.polyfit(FFup,chiup,degree) #yfit in ron's code
        chifit_v = np.polyval(chifit,FFup)
        chi = chiup - chifit_v
        
        [toss, unique_ind] = np.unique(FFinv, return_index = True)
        unique_ind = np.sort(unique_ind).astype(int)
        FF = (np.array(FFup))[unique_ind]
        chi = (np.array(chiup))[unique_ind]
        
        chispl_fn = interpolate.interp1d(1/FF, chi)  #np.interp(FFinvspl,1/FF,chi)
        chispl = chispl_fn(FFinvspl);
        n_datapoints = len(chispl) #Datapoints in ron's code
        hwind = np.hanning(n_datapoints)
        chisplW = chispl*hwind; #ysplW in ron's code
        
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
        fI = np.where(f> self.inputFrequencyRange[0] & f < self.inputFrequencyRange[1])
        f = f[fI]
        
        self.FFTCalcObj.ffInv = 1/FF
        self.FFTCalcObj.chi = chi
        self.FFTCalcObj.chiSpline = chispl
        self.FFTCalcObj.ffInvSpline = FFinvspl
        self.FFTCalcObj.frequencies = f
        self.FFTCalcObj.FFT = FFTval[fI]
        self.FFTCalcObj.phase = FFTphase[fI]
        
        print('Finished calculating the FFT\n')
        
    def func_FindPeaks(self,minFreq = 150): 
        FFT = self.FFTCalcObj.FFT
        f = self.FFTCalcObj.frequencies
        rms = np.std(FFT)
        peakLoc = signal.find_peaks(FFT,height = rms)[0]
        peakVal = FFT[peakLoc]
        peakFreq = f[peakLoc]
        peakValReduced = []
        peakFreqReduced = []
        
        for ii in range(len(peakVal)):
            nearbyPeakLocs = np.where(abs(peakFreq[ii] - peakFreq) < 150)
            nearbyPeakVals = peakVal[nearbyPeakLocs]
            nearbyPeakFreqs = peakFreq[nearbyPeakLocs]
            reduced_peak_loc = np.where(peakVal == max(nearbyPeakVals))[0]
            if len(peakFreqReduced) != 0 and min(nearbyPeakFreqs) > minFreq:
                if max(nearbyPeakVals) != peakValReduced[-1]:
                    peakValReduced.append(max(nearbyPeakVals))
                    peakFreqReduced.append(peakFreq[reduced_peak_loc][0])
            elif len(peakFreqReduced) == 0 and min(nearbyPeakFreqs) > minFreq:
                    peakValReduced.append(max(nearbyPeakVals))
                    peakFreqReduced.append(peakFreq[reduced_peak_loc][0])
                    
        return [peakFreq,peakVal]
    
    def func_MassCalc(self):
        
        if self.FFTCalcObj.frequencies == []:
            print('Calculating FFT first\n')
            self.FFTCalcObj(self)
        [peakFreqs,peakAmps] = self.func_FindPeaks()
        f = self.FFTCalcObj.frequencies
        FFTvs = self.FFTCalcObj.FFT
        tempVec = []
        if len(self.MassCalcObj.tempVec) == 0:
            print("A vector of temperatures is neeeded for fitting.\n")
            print("Assign a vector of temperatures to MassCalcObject.tempVec.")
            exit()
        elif len(self.MassCalcObj.tempVec) != 0:
            tempVec = self.MassCalcObj.tempVec
            
        
        A = max(FFTvs)
        maxI = FFTvs.argmax()
        temp = self.inputTemperature
        Brange = self.FFTCalcObj.FFrange
        MF = f[maxI]
        maxFreq = np.mean(MF)
        phaseV = self.FFTCalcObj.phase
        phase = phaseV[maxI]
        AoT = A/temp
        
        def fit_func(a,b,x):
            return  a/np.sinh(b*x)
        Bm= 1/(0.5*(1/Brange[0]) + 1/Brange[-1])
        fitresult = optimize.least_squares(fit_func)
            
        return print('MassCalc')
        

    
