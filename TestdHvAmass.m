clear
close all
clc
%% Here the files to be loaded are chosen.
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\031919_PdGa\031919\';
files = dir(strcat(loc,'*.ASC'));
filesFFT = dir(strcat(loc,'FT*'));
temps = [1.52 .62 5.71 3.96 8 10 12 14 18 3.27 2.24 .96];%Temperatures
% corresponding to files to be loaded
I = ismember(temps,temps);%I is the index that will load all files 
filesa = files(I);
filesFFTa = filesFFT(I);
T1 = temps(1:2:end);
I1 = ismember(temps,T1);
files1 = files(I1);
filesFFT1 = filesFFT(I1);
T2 = temps(2:2:end);
I2 = ismember(temps,T2);
files2 = files(I2);
filesFFT2 = filesFFT(I2);

%% Here the raw data is loaded into dHvA class
obj = dHvA(loc,temps,filesa,filesFFTa);

%% Here the 1/H windows that are selected and fourier transformed are defined
endFields = 35;%the endFields are the maximum field values of each window
iFFspan = abs((1/35-1/15));%Here the width of the 1/H window is defined 

%% Here the Fourier transform is performed  
dHvA.FFTload(obj,endFields,iFFspan,[5 4500]);

%% 
peakRange = [490 708; 708 850;  2600 2900; 3546 3880]; 
dHvA.massLoad(obj,peakRange);

%% Plot AoT vs. temp from my MATLAB fft
for jj = 1:length(obj.FFT.range)
    figure
    for ii = 1:length(obj.mass.range(1).upPeak)
        
        subplot(2,2,ii)
        plot(obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoT,'*',...
            obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoTfitV)
%         maxFreq = obj.mass.range.upPeak(kk).maxFreq;
        title('MATLAB All Temps')
        xlabel('Temperature (K)')
        ylabel('Amplitude/T (arb. units)')
        
    end
    suptitle('Fig 3.a: MATLAB ft, Fourier Window 15 to 35 T')
end
%% AoT vs temp plot from LANL fft
figure
for kk = 1:4
    subplot(2,2,kk)
    plot(obj.mass.raw(kk).temp,obj.mass.raw(kk).AoT,'*',...
        obj.mass.raw(kk).temp,obj.mass.raw(kk).AoTfitV)
%     maxFreq = obj.mass.range.upPeak(kk).maxFreq;
    title('LANL All Temps')
    xlabel('Temperature (K)')
    ylabel('Amplitude/T (arb. units)')
    
end

suptitle('Fig 3.b: LANL ft, Fouriet Window 15 to 35 T')

%% Here the raw data is loaded into dHvA class
obj1 = dHvA(loc,files1,filesFFT1,T1);

%% Here the 1/H windows that are selected and fourier transformed are defined
endFields = 35;%the endFields are the maximum field values of each window
iFFspan = abs((1/35-1/15));%Here the width of the 1/H window is defined 

%% Here the Fourier transform is performed  
dHvA.FFTload(obj1,endFields,iFFspan,[5 4500]);

%% Calculate effective masses 
peakRange = [490 708; 708 850;  2600 2900; 3546 3880]; 
dHvA.massLoad(obj1,peakRange);
% Plot AoT vs. temp from my MATLAB fft
for jj = 1:length(obj1.FFT.range)
    figure
    for ii = 1:length(obj1.mass.range(1).upPeak)
        
        subplot(2,2,ii)
        plot(obj1.mass.range(jj).upPeak(ii).temp,obj1.mass.range(jj).upPeak(ii).AoT,'*',...
            obj1.mass.range(jj).upPeak(ii).temp,obj1.mass.range(jj).upPeak(ii).AoTfitV)
%         maxFreq = obj.mass.range.upPeak(kk).maxFreq;
        title('MATLAB obj1')
        xlabel('Temperature (K)')
        ylabel('Amplitude/T (arb. units)')
        
    end
    suptitle('MATLAB ft, Fourier Window 15 to 35 T')
end
%% AoT vs temp plot from LANL fft
figure
for kk = 1:4
    subplot(2,2,kk)
    plot(obj1.mass.raw(kk).temp,obj1.mass.raw(kk).AoT,'*',...
        obj1.mass.raw(kk).temp,obj1.mass.raw(kk).AoTfitV)
%     maxFreq = obj.mass.range.upPeak(kk).maxFreq;
    title('LANL obj1')
    xlabel('Temperature (K)')
    ylabel('Amplitude/T (arb. units)')
    
end

suptitle('LANL ft, Fouriet Window 15 to 35 T')

%% Here the raw data is loaded into dHvA class
obj2 = dHvA(loc,files1,filesFFT1,T2);

%% Here the 1/H windows that are selected and fourier transformed are defined
endFields = 35;%the endFields are the maximum field values of each window
iFFspan = abs((1/35-1/15));%Here the width of the 1/H window is defined 

%% Here the Fourier transform is performed  
dHvA.FFTload(obj2,endFields,iFFspan,[5 4500]);

%% Calculate effective masses 
peakRange = [490 708; 708 850;  2600 2900; 3546 3880]; 
dHvA.massLoad(obj2,peakRange);
%% Plot AoT vs. temp from my MATLAB fft
for jj = 1:length(obj1.FFT.range)
    figure
    for ii = 1:length(obj1.mass.range(1).upPeak)
        
        subplot(2,2,ii)
        plot(obj2.mass.range(jj).upPeak(ii).temp,obj2.mass.range(jj).upPeak(ii).AoT,'*',...
            obj2.mass.range(jj).upPeak(ii).temp,obj2.mass.range(jj).upPeak(ii).AoTfitV)
%         maxFreq = obj.mass.range.upPeak(kk).maxFreq;
        title('MATLAB obj1')
        xlabel('Temperature (K)')
        ylabel('Amplitude/T (arb. units)')
        
    end
    suptitle('MATLAB ft, Fourier Window 15 to 35 T')
end
%% AoT vs temp plot from LANL fft
figure
for kk = 1:4
    subplot(2,2,kk)
    plot(obj2.mass.raw(kk).temp,obj2.mass.raw(kk).AoT,'*',...
        obj2.mass.raw(kk).temp,obj2.mass.raw(kk).AoTfitV)
%     maxFreq = obj.mass.range.upPeak(kk).maxFreq;
    title('LANL obj2')
    xlabel('Temperature (K)')
    ylabel('Amplitude/T (arb. units)')
    
end

suptitle('LANL ft, Fouriet Window 15 to 35 T')

%% Plot masses vs frequency
figure

for jj = 1:length(obj.FFT.range)
%     jj
    for ii = 1:length(obj.mass.range(1).upPeak)
        massPlot(ii) = obj.mass.range(jj).upPeak(ii).m;
        freqPlot(ii) = mean(obj.mass.range(jj).upPeak(ii).maxFreq);
        rawMassPlot(ii) = obj.mass.raw(ii).m;
        rawFreqPlot(ii) = mean(obj.mass.raw(ii).maxFreq);
    end
    
    for ii = 1:length(obj1.mass.range(1).upPeak)
        massPlot1(ii) = obj1.mass.range(jj).upPeak(ii).m;
        freqPlot1(ii) = mean(obj1.mass.range(jj).upPeak(ii).maxFreq);
        rawMassPlot1(ii) = obj1.mass.raw(ii).m;
        rawFreqPlot1(ii) = mean(obj1.mass.raw(ii).maxFreq);
    end
    
    for ii = 1:length(obj2.mass.range(1).upPeak)
        massPlot2(ii) = obj2.mass.range(jj).upPeak(ii).m;
        freqPlot2(ii) = mean(obj2.mass.range(jj).upPeak(ii).maxFreq);
        rawMassPlot2(ii) = obj2.mass.raw(ii).m;
        rawFreqPlot2(ii) = mean(obj2.mass.raw(ii).maxFreq);
    end
    leg{jj} = num2str(obj.FFT.range(jj).upTemp(1).range);
%     errorbar(freqPlot,massPlot,errorPlot,'o')
    plot(freqPlot,massPlot,'o',rawFreqPlot,rawMassPlot,'*')
    pause(1)
    hold on
    plot(freqPlot1,massPlot1,'o',rawFreqPlot1,rawMassPlot1,'*')
    pause(1)
    plot(freqPlot2,massPlot2,'o',rawFreqPlot2,rawMassPlot2,'*')

end
% errorbar(rawFreqPlot,rawMassPlot,rawErrorPlot,'*')
title('Fourier Window 15 to 35')
xlabel('Effective mass, m*')
ylabel('F(T)')
legend('m* from MATLAB ft','m* from LANL ft','m* from MATLAB ft1',...
    'm* from LANL ft1','m* from MATLAB ft2','m* from LANL ft2')
