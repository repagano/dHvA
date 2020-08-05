import numpy as np
from scipy import interpolate,signal, optimize
import pandas as pd
import os

#%%
def dHvALoad(fn,FFmin = 1.5):
    data = pd.read_csv(fn, delimiter = '\t')
    Bdot = np.array(data.Bdot) #dH/dt
    Mdot = np.array(data.Mag1) #dM/dt (Mag1 in ron's code)
    FF = np.array(data.Field_fixed)
    chi = Mdot/Bdot #(xi in ron's code)
    FFmax_ind = FF.argmax()
    FF[FF < FFmin] = 0
    FFup = FF[0:FFmax_ind]
    chiup = [];
    FFup_return = []
    for ii in range(len(FFup)):
        if FFup[ii] != 0:
            chiup.append(chi[ii])
            FFup_return.append(FFup[ii])
    return chiup, FFup_return
    

def FFTCalc(fn,FFmin = 1.5, FFrange = [10,50], freqRange = [400,4000]):
    [chi,FF] = dHvALoad(fn = fn, FFmin = FFmin)
    FFupi = np.array(FF)
    chiupi = np.array(chi)
    FFupi[FFupi < FFrange[0]] = 0
    FFupi[FFupi > FFrange[1]] = 0
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
    FFTPhases = FFTphase[0:int(n_datapoints/2 + 1)]
    FFTAmplitudes = abs(FFTcomplex[0:int(n_datapoints/2+1)])
    
    length = abs(ffinv_fft[-1] - ffinv_fft[1])
    fs = n_datapoints/length #sampling frequency
    f = fs/2 * np.linspace(0,1,n_datapoints/2 + 1); #frequency
    fI = np.where((f> freqRange[0]) & (f < freqRange[1]))
    frequencies = f[fI]
    
    return frequencies, FFTAmplitudes[fI], FFTPhases[fI]

        
def findPeaks(fn = '',FFT = [],f = [],freqRange = [400,4000], FFrange = [10,50],FFmin = 1.5):
    if fn != '':
        if FFT != []:
            print('Please only provide a filename OR FFT and frequency data, not both.\n')
            exit(keep_kernel = True);
        elif FFT == []:
            [f,FFT,toss] = FFTCalc(fn,FFmin,FFrange,freqRange)
            
    elif fn == '':
        if len(FFT) != len(f):
            print('Length of FFT and f must match.\n')
            exit(keep_kernel = True)
    
    rms = np.std(FFT)
    peakLoc = signal.find_peaks(FFT,height = rms)[0]
    peakVal = FFT[peakLoc]
    peakFreq = f[peakLoc]
    peakValReduced = []
    peakFreqReduced = []
    minFreq = freqRange[0];
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
                
    return peakFreqReduced, peakValReduced


class MassCalc():
    fnVec= []
    tempVec = []
    freqRangeVec = []
    FFRangeVec = []
    FFmin = 1.5
    
    _Processed = False
    
    _col1 = 'Temperature'
    _col2 = 'Filepath'
    _col3 = 'Frequencies'
    _col4 = 'FFT_Amplitudes'
    _col5 = 'FFT_Phases'
    _col6 = 'Peak_Locations'
    _col7 = 'Peak_Amplitudes'
    _col8 = 'Peak_Phases'
    _col9 = 'Amplitude_rms'
    _dataFrameHeaderTemp = [_col1,_col2,_col3,_col4,_col5,_col6,_col7,_col8,_col9]
    dataFrameTemp = pd.DataFrame(columns = _dataFrameHeaderTemp)
    
    _col1n = 'Average_Peak_Frequency'
    _col2n = 'AoT'
    _col3n = 'Mass'
    _col4n = 'Mass_Calc_Error'
    _col5n = 'AOT_Parameters'
    _dataFrameHeaderPeak = [_col1n,_col2n,_col5n,_col3n,_col4n]
    dataFramePeak = pd.DataFrame(columns = _dataFrameHeaderPeak)
    
    def _checkProc(self):
        try:
            self.dataFrameTemp.iloc[0]
            self._Processed = True
        except:
            self._Processed = False
        
    
    def _dataPrep(self):
        #Putting the temperatures and corresponding fn's in ascending order:
        [newTempVec,newTempVecLocs] = np.unique(self.tempVec,return_index = True)
        self.tempVec = newTempVec
        self.fnVec = np.array(self.fnVec)[newTempVecLocs]
        
        #Updating table of info for each temperature:
        for ii in range(len(self.fnVec)):
            [frequencies, FFTAmplitudes, FFTPhases] = FFTCalc(self.fnVec[ii],self.FFmin,self.FFRangeVec[ii],self.freqRangeVec[ii]);
            col1 = self.tempVec[ii]
            col2 = self.fnVec[ii];
            col3 = frequencies
            col4 = FFTAmplitudes
            col5 = FFTPhases
            col6 = [0]
            col7 = [0]
            col8 = [0]
            col9 = [0]
            dat = {self._col1:col1, self._col2:col2, self._col3:col3, self._col4:col4, self._col5:col5,self._col6:col6, self._col7:col7, self._col8:col8, self._col9:col9}
            self.dataFrameTemp = self.dataFrameTemp.append(dat,ignore_index = True)
       
    def _ComparePeaks(self):
        peakFreqs = []
        peakAmps = []
        significantPeaks = np.empty((0,0))
        for ii in range(len(self.tempVec)):
            [freqs, amps] = findPeaks(fn = self.fnVec[ii],freqRange=self.freqRangeVec[ii],FFrange=self.FFRangeVec[ii],FFmin=self.FFmin)
            peakFreqs.append(freqs); 
            peakAmps.append(amps)
        
        for ii in range(len(self.tempVec)):
            currentFreqs = np.array(peakFreqs[ii])
            for jj in range(len(self.tempVec)):
                testFreq = np.array(peakFreqs[ii])
                if ii != jj:
                    for kk in range(len(currentFreqs)):
                        simFreqLoc = np.where(abs(currentFreqs[kk]-testFreq) <= 10)[0]
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

        
    def _findCorrespondingAmp(self,ind,commonPeakFreqs = []):
        if commonPeakFreqs == []:
            commonPeakFreqs = self._ComparePeaks()
        FFTObj = self.dataFrameTemp.iloc[ind]
        FFT = np.array(FFTObj.FFT_Amplitudes)
        f = np.array(FFTObj.Frequencies)
        phases = np.array(FFTObj.FFT_Phases)
        peakAmp = []
        rmsAmp = []
        peakFreq = []
        peakPhase = []
        for freqPeak in commonPeakFreqs:
            freqWindow = 0.02*freqPeak
            fWindowMin = freqPeak - freqWindow
            fWindowMax = freqPeak + freqWindow
            searchLoc = np.where((f > fWindowMin) & (f < fWindowMax))[0]
            
            ampVec = FFT[searchLoc]
            maxAmp = ampVec.max()
            maxInd = np.where(FFT == maxAmp)
            maxFreq = f[maxInd]
            maxPhase = phases[maxInd]
            
            newWindow = [maxFreq - maxFreq*0.075, maxFreq + maxFreq*0.075]
            newWindowInd = np.where( (f > newWindow[0]) & (f < newWindow[1]))
            newWindowAmp = FFT[newWindowInd]
            rms = np.std(newWindowAmp)
            
            peakAmp.append(maxAmp)
            peakFreq.append(maxFreq[0])
            rmsAmp.append(rms)
            peakPhase.append(maxPhase[0])
        self.dataFrameTemp.Peak_Locations[ind] = peakFreq
        self.dataFrameTemp.Peak_Amplitudes[ind] = peakAmp
        self.dataFrameTemp.Peak_Phases[ind] = peakPhase
        self.dataFrameTemp.Amplitude_rms[ind] = rmsAmp

        
    def func_MassCalc(self):
        self._checkProc()
        if self._Processed == False:
            self._dataPrep()
        commonPeakFreqs = self._ComparePeaks()
        tempVec = np.array(self.tempVec)
        tempLength = len(tempVec)
        freqLength = len(commonPeakFreqs)
        frequencyMat = np.empty((tempLength,freqLength))#horizontal axis: frequency, vertical axis: temperature
        AoTMat = np.empty((tempLength,freqLength))
        rmsMat = np.empty((tempLength,freqLength))
        BmVec = []
        
        for ii in range(len(tempVec)):
            self._findCorrespondingAmp(ii,commonPeakFreqs)
            peakFreq = self.dataFrameTemp.Peak_Locations[ii]
            peakAmp =  self.dataFrameTemp.Peak_Amplitudes[ii]
            rmsAmp = self.dataFrameTemp.Amplitude_rms[ii]
            frequencyMat[ii] = peakFreq
            AoTMat[ii] = peakAmp/tempVec[ii]
            rmsMat[ii] = rmsAmp
            FFrange = self.FFRangeVec[ii]
            BmVec.append(1/(0.5*(1/FFrange[0] + 1/FFrange[-1])))
            
        Bm = np.unique(BmVec)[0]
        p0 = [0,0.1]
        
        def AoTFitType(x,a,b):
            return a/np.sinh(b*x)

        for ii in range(len(AoTMat.transpose())):
            dataToFit = AoTMat[:,ii] 
            AoTfit_params, params_covar = optimize.curve_fit(AoTFitType,tempVec,dataToFit,p0 = p0,sigma = rmsMat[:,ii])
            col1nval = np.average(frequencyMat[:,ii])
            col2nval = dataToFit
            col3nval = AoTfit_params[1]*Bm/14.69
            col4nval = np.sqrt(np.diag(params_covar))[1]
            col5nval = AoTfit_params
            datn = {self._col1n:col1nval,self._col2n:col2nval,self._col3n:col3nval,self._col4n:col4nval,self._col5n:col5nval}
            self.dataFramePeak = self.dataFramePeak.append(datn,ignore_index = True)
