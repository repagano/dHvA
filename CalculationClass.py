import numpy as np
#import scipy as sp
from scipy import interpolate,signal, optimize
#import matplotlib.pyplot as plt
import pandas
#import csv
#import os

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
    

class dataProcess(): 
    inputFrequencyRange = []
    fn = ''
    inputTemperature = []
    dHvALoadObj = ClassdHvALoad()
    FFTCalcObj = ClassFFTCalc()
    #MassCalcObj = ClassMassCalc()
    
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
        self.dHvALoadObj.FFup = np.array(FFup_return)
        #print('Finished loading dHvA data\n')
        
    def func_FFTCalc(self):
        if self.dHvALoadObj.FFup != [] and self.dHvALoadObj.chiup != []:
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
        nonzeroIndices = np.where(FFupi != 0)[0]
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
        fI = np.where((f> self.inputFrequencyRange[0]) & (f < self.inputFrequencyRange[1]))
        f = f[fI]
        
        self.FFTCalcObj.ffInv = 1/FF
        self.FFTCalcObj.chi = chi
        self.FFTCalcObj.chiSpline = chispl
        self.FFTCalcObj.ffInvSpline = FFinvspl
        self.FFTCalcObj.frequencies = f
        self.FFTCalcObj.FFT = FFTval[fI]
        self.FFTCalcObj.phase = FFTphase[fI]
        
        print('Finished calculating the FFT\n')
        
#    def func_FindPeaks(self,minFreq = 150): 
#        FFT = self.FFTCalcObj.FFT
#        f = self.FFTCalcObj.frequencies
#        rms = np.std(FFT)
#        peakLoc = signal.find_peaks(FFT,height = rms)[0]
#        peakVal = FFT[peakLoc]
#        peakFreq = f[peakLoc]
#        peakValReduced = []
#        peakFreqReduced = []
#        
#        for ii in range(len(peakVal)):
#            nearbyPeakLocs = np.where(abs(peakFreq[ii] - peakFreq) < 150)
#            nearbyPeakVals = peakVal[nearbyPeakLocs]
#            nearbyPeakFreqs = peakFreq[nearbyPeakLocs]
#            reduced_peak_loc = np.where(peakVal == max(nearbyPeakVals))[0]
#            if len(peakFreqReduced) != 0 and min(nearbyPeakFreqs) > minFreq:
#                if max(nearbyPeakVals) != peakValReduced[-1]:
#                    peakValReduced.append(max(nearbyPeakVals))
#                    peakFreqReduced.append(peakFreq[reduced_peak_loc][0])
#            elif len(peakFreqReduced) == 0 and min(nearbyPeakFreqs) > minFreq:
#                    peakValReduced.append(max(nearbyPeakVals))
#                    peakFreqReduced.append(peakFreq[reduced_peak_loc][0])
#                    
#        return [peakFreq,peakVal]
       
#    def func_MassCalc(self):
#        
#        if self.FFTCalcObj.frequencies == []:
#            print('Calculating FFT first\n')
#            self.FFTCalcObj(self)
#        [peakFreqs,peakAmps] = self.func_FindPeaks()
#        f = self.FFTCalcObj.frequencies
#        FFTvs = self.FFTCalcObj.FFT
#        tempVec = []
#        if len(self.MassCalcObj.tempVec) == 0:
#            print("A vector of temperatures is neeeded for fitting.\n")
#            print("Assign a vector of temperatures to MassCalcObject.tempVec.")
#            exit()
#        elif len(self.MassCalcObj.tempVec) != 0:
#            tempVec = np.array(self.MassCalcObj.tempVec)
#            
#        
#        A = max(FFTvs)
#        maxI = FFTvs.argmax()
#        temp = self.inputTemperature
#        Brange = self.FFTCalcObj.FFrange
#        MF = f[maxI]
#        maxFreq = np.mean(MF)
#        phaseV = self.FFTCalcObj.phase
#        phase = phaseV[maxI]
#        AoT = A/temp
#        
#        Bm= 1/(0.5*(1/Brange[0]) + 1/Brange[-1])
#        fitresult = optimize.least_squares(lambda coeffs:coeffs[0]/np.sinh(coeffs[1]*tempVec) )
#            
#        return print('MassCalc')
        
class MassCalc():
    fnVec= []
    tempVec = []
    freqRangeVec = []
    massCalcObj = ClassMassCalc();
    dataProcObjVec = []
#    dHvALoadObj = dataProcess.dHvALoadObj();
#    FFTCalcObj = dataProcess.FFTCalcObj();
    _dataProcVecNonzero = False
    
    def _checkProcVec(self):
        if len(self.dataProcObjVec) == 0:
            self._dataProcVecNonzero = False;
        else:
            self._dataProcVecNonzero = True
        return self._dataProcVecNonzero 
        
    
    def func_dataPrep(self):
        dataProcBool = self._checkProcVec();
        if dataProcBool == False:
            if len(self.fnVec) != len(self.tempVec) or len(self.fnVec) != len(self.freqRangeVec) :
                print("fnVec and tempVec must be the same length")
                exit();
            for ii in range(len(self.fnVec)):
                dataProcObj = dataProcess();
                dataProcObj.fn = self.fnVec[ii]
                dataProcObj.inputTemperature = self.tempVec[ii]
                dataProcObj.inputFrequencyRange = self.freqRangeVec[ii]
                dataProcObj.func_dHvALoad()
                dataProcObj.func_FFTCalc()
                self.dataProcObjVec.append(dataProcObj)
        
    def _findPeaks(self,index):
        FFTObj = self.dataProcObjVec[index].FFTCalcObj;
        FFT = FFTObj.FFT
        f = FFTObj.frequencies
        rms = np.std(FFT)
        peakLoc = signal.find_peaks(FFT,height = rms)[0]
        peakVal = FFT[peakLoc]
        peakFreq = f[peakLoc]
        peakValReduced = []
        peakFreqReduced = []
        minFreq = self.dataProcObjVec[index].inputFrequencyRange[0];
        for ii in range(len(peakVal)):
            nearbyPeakLocs = np.where(abs(peakFreq[ii] - peakFreq) < 50)
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
    
    def _ComparePeaks(self):
        peakFreqs = []
        peakAmps = []
        significantPeaks = np.empty((0,0))
        for ii in range(len(self.dataProcObjVec)):
            [freqs,amps] = self._findPeaks(self,ii)
            peakFreqs.append(freqs)
            peakAmps.append(amps)
        
        for ii in range(len(self.dataProcObjVec)):
            currentFreqs = np.array(peakFreqs[ii])
            for jj in range(len(self.dataProcObjVec)):
                testFreq = np.array(peakFreqs[ii])
                if ii != jj:
                    for kk in range(len(currentFreqs)):
                        simFreqLoc = np.where(abs(currentFreqs[kk]-testFreq) <= 10)
                        if simFreqLoc.size != 0:
                            significantPeaks = np.append(significantPeaks, testFreq[simFreqLoc])
        significantPeaks = np.unique(significantPeaks)
        
        refinedFrequencies = np.empty((0,0))
        for ii in range(len(significantPeaks)):
            currentFreq = significantPeaks[ii]
            simFreqLocs = np.where(abs(currentFreq - significantPeaks)<9)[0]
            if simFreqLocs.size > 5:
                refinedFrequencies = np.append(refinedFrequencies,significantPeaks[simFreqLocs])
        
        refinedFrequencies = np.unique(refinedFrequencies)
        refinedFrequencies = list(refinedFrequencies)
        
        ind = 0
        currentFreq = refinedFrequencies[0]
        while currentFreq < refinedFrequencies[-1]:
            currentFreq = refinedFrequencies[ind]
            oom = np.floor(np.log10(currentFreq))
            ind+=1
            if oom <= 2.:
                fact = 0.1
            else:
                fact = 0.02
            diff = abs(np.array(refinedFrequencies)-currentFreq)
            badLocs = np.where((diff < fact*currentFreq) & (refinedFrequencies != currentFreq))
            badFreqs = np.array(refinedFrequencies)[badLocs]
            if len(badFreqs) != 0:
                for badFreq in badFreqs:
                    refinedFrequencies.remove(badFreq)
                ind = 0
        return np.unique(refinedFrequencies)
#        finalFrequencies = []
#        for ii in range(len(refinedFrequencies)):
#            currentFreq = refinedFrequencies[ii]
#            if ii != 0:
#                if currentFreq - refinedFrequencies[ii-1] < 20:
#                    finalFrequencies.append((currentFreq + refinedFrequencies[ii-1])/2)
#                else:
#                    finalFrequencies.append(currentFreq)
#            if ii == 0:
#                finalFrequencies.append(currentFreq)
#        return np.unique(np.array(finalFrequencies));
        
    def _findCorresponingAmp(self,ind):
        commonPeakFreqs = self._ComparePeaks()
        FFTObj = self.dataProcObjVec[ind].FFTCalcObj
        FFT = np.array(FFTObj.FFT)
        f = np.array(FFTObj.frequencies)
        peakAmp = []
        rmsAmp = []
        peakFreq = []
        #want to find where the max value occurs and then do a window around IT
        for freqPeak in commonPeakFreqs:
            #oom = np.floor(np.log10(freqPeak))
            freqWindow = 0.02*freqPeak
            fWindowMin = freqPeak - freqWindow
            fWindowMax = freqPeak + freqWindow
            searchLoc = np.where((f > fWindowMin) & (f < fWindowMax))[0]
            ampVec = FFT[searchLoc]
            maxAmp = ampVec.max()
            maxInd = np.where(FFT = maxAmp)
            maxFreq = f[maxInd]
            rms = np.std(ampVec)
            peakAmp.append(maxAmp)
            peakFreq.append(maxFreq)
            rmsAmp.append(rms)
        return peakFreq, peakAmp, rmsAmp
                
            
        
            
    

        
    def func_MassCalc(self):
        if self._checkProcVec == False:
            self.func_dataPrep()
        for dataProcObj in self.dataProcObjVec:
            FFTObj = dataProcObj.FFTCalcObj;
            
        
        
                
    
            
 
folder_path = r'C:\Users\iveli\OneDrive\Documents\Ron\031919'
file_paths = []

for file in os.listdir(folder_path):
    if file.endswith('.ASC'):
        file_paths.append(folder_path+'\\'+str(file))
        
temps = [1.52, 0.62, 5.71, 3.96, 8, 10, 12, 14, 18, 3.27, 2.24, 0.96];    

#testObj = dataProcess()
#testObj.fn = file_paths[0]
#testObj.inputTemperature = temps[0]
#testObj.inputFrequencyRange = [10,8000]
#
#testObj.func_dHvALoad()
#testObj.func_FFTCalc()

masstestObj = MassCalc()
masstestObj.fnVec = file_paths
masstestObj.tempVec = temps
frang= [400,4000]
masstestObj.freqRangeVec = [frang]*len(file_paths)
masstestObj.func_dataPrep()