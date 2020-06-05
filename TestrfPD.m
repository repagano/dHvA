clear
close all
clc
%% Here the files to be loaded are chosen.
addpath C:\Users\repag\Documents\MATLAB\DiTusa
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\CoSi\';
temps = [4.03 3.98 3.31 2.75 2.1 1.57 0.621 0.56 8 12 10 6];%Temperatures
% corresponding to files to be loaded
files = dir(loc);
% filesFFT = dir(strcat(loc,'FT*'));
% here set the variable temps equal to the data temperatures corresponding 
% to the files you wish to load in the expression find(temps == temps) by 
% setting temps == temps, all files are loaded.
I = ismember(temps,temps(2:12));
% I = logical([0 I(2:end)]);
temps = temps(I);
files = files(logical([0, 0 , I]));
% filesFFT = filesFFT(I);

%% Here the raw data is loaded into dHvA class
obj = dHvA(loc,temps,files);

%% Here the 1/H windows that are selected and fourier transformed are defined
endFields = [20:60];%the endFields are the maximum field values of each window
iFFspan = abs((1/22-1/60));%Here the width of the 1/H window is defined 

%% Here the Fourier transform is performed  
dHvA.FFTload(obj,[22,60],iFFspan,[5 1000]);

%% Here the effective masses of the different peaks are calculated 
% peakRange is an matrix whose rows define the range aro und each peak 
% identified on the Fourier spectrum. The columns correspond to different
% peaks. The range should be large enough to allow for variation in the position 
% of the peak as the different ranges and temperatures are scanned over, but 
% not be so large that another peak is selected. 
peakRange = [158 565; 565 680]; 
% dHvA.massLoad runs the method that calculates the effective masses.
dHvA.massLoad(obj,peakRange);

%% plot dM/dH vs 1/H. Fig 1 in document created using this section with
% endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
temp = 1:length(obj.raw);
for jj = length(obj.FFT.range)
    for ii = temp
        
        figure
%         subplot(2,1,1)
        plot(obj.FFT.range(jj).upTemp(ii).x,obj.FFT.range(jj).upTemp(ii).y,'Linewidth',1.2)
        hold on
        plot(obj.FFT.range(jj).upTemp(ii).xspl,obj.FFT.range(jj).upTemp(ii).yspl,'Linewidth',.4)
        xlabel('1/H')
%         ylabel('dM/dH')
%         legend('dM/dH vs 1/H','Spline fit')
%         subplot(2,1,2)
%         plot(1./obj.FFT.range(jj).upTemp(ii).x,obj.FFT.range(jj).upTemp(ii).y,'Linewidth',1)
%         hold on
%         plot(1./obj.FFT.range(jj).upTemp(ii).xspl,obj.FFT.range(jj).upTemp(ii).yspl,'Linewidth',.4)
%         xlabel('H')
%         ylabel('dM/dH')
%         legend('dM/dH vs 1/H','Spline fit')
    end      
end
suptitl1 = strcat(num2str(obj.raw(ii).temp),'K');
suptitle(strcat('Fig 1: T=',suptitl1))

%% Fourier transform from MATLAB.
temps = 1:length(obj.raw);
pos = 1;
figure
for ii = temps

%     subplot(2,2,pos)
%     plot(obj.raw(ii).f,obj.raw(ii).FFT)
%     hold on
%     figure
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).FFT)
    hold on
    leg{ii} = strcat(num2str(obj.raw(ii).temp),'K');
    title(leg{ii})
%     ylim([0 11])
%     legend('LANL fourier transform','MATLAB fourier transform')
%     xlabel('dHvA frequency, F (T)')
%     ylabel('Amplitude (arb. units)')
    pos = pos+1;
%     for jj = 1:4
%         plot(obj.mass.raw(jj).maxFreq(ii),obj.mass.raw(jj).A(ii),'*')
%     end
end
title('fourier window 22 to 60 T')
legend(leg)
 
%% Plot max peaks vs. temps
figure
for ii = 1:length(obj.mass.range.upPeak)
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).A,'-*')
    hold on
    leg{ii} = strcat('f = ', num2str(mean(obj.mass.range.upPeak(ii).maxFreq)));
end
xlabel('T (K)')
ylabel('Amplitude (arb. units)')
legend(leg) 
%% Plot masses vs frequency
figure

for jj = 1:length(obj.FFT.range)
%     jj
    for ii = 1:length(obj.mass.range(1).upPeak)
        massPlot(ii) = obj.mass.range(jj).upPeak(ii).m;
        freqPlot(ii) = mean(obj.mass.range(jj).upPeak(ii).maxFreq);
%         errorPlot(ii) = mean(obj.mass.range(jj).upPeak(ii).AoTrms);
%         rawMassPlot(ii) = obj.mass.raw(ii).m;
%         rawFreqPlot(ii) = mean(obj.mass.raw(ii).maxFreq);
%         rawErrorPlot(ii) = mean(obj.mass.raw(ii).AoTrms);
    end
    leg{jj} = num2str(obj.FFT.range(jj).upTemp(1).range);
%     errorbar(freqPlot,massPlot,errorPlot,'o')
    plot(freqPlot,massPlot,'o')%rawFreqPlot,rawMassPlot,'*')
    hold on

end
% errorbar(rawFreqPlot,rawMassPlot,rawErrorPlot,'*')
title('Fourier Window 15 to 35')
xlabel('Effective mass, m*')
ylabel('F(T)')
legend(leg)%'m* from MATLAB ft','m* from LANL ft')

%% plot Fourier transform for all temps 
cnt = [1 2];% 3 4 5 6 7 10];
sp = 1;
figure
for jj = 1:length(obj.FFT.range)
%     subplot(2,3,sp)
    for ii = 1:length(obj.raw)
%         figure(jj)
        leg{ii} = [strcat(num2str(obj.FFT.range(jj).upTemp(ii).temp), 'K')];
%         subplot(2,4,sp)
        plot(obj.FFT.range(jj).upTemp(ii).f,obj.FFT.range(jj).upTemp(ii).FFT)
        hold on
        title(num2str(obj.FFT.range(jj).upTemp(ii).range))
    
    end

sp = sp+1;
titl1 = strcat('CF = ',num2str(round(mean(obj.FFT.range(jj).upTemp(1).range))));
titl = strcat(titl1,'T');
title(titl)
xlabel('frequency (T)')
ylabel('Amplitude (arb. units)')

end
legend(leg)

%% plot fourier transform for temps individually 
cnt = [1 2];% 2 3 4 5 6 7 10];
sp = 1;
figure
for ii = 1:length(obj.FFT.range)%cnt
    
%     subplot(2,3,sp)    
    plot(obj.FFT.range(ii).upTemp(2).f,obj.FFT.range(ii).upTemp(2).FFT,'r')
    hold on
    plot(obj.FFT.range(ii).upTemp(3).f,obj.FFT.range(ii).upTemp(3).FFT,'b')
    plot(obj.FFT.range(ii).upTemp(7).f,obj.FFT.range(ii).upTemp(7).FFT,'Color',[0 .7 .3])
    titl1 = strcat('CF = ',num2str(round(mean(obj.FFT.range(ii).upTemp(1).range))));
    titl = strcat(titl1,'T');
    title(titl)
    xlabel('frequency (T)')
    ylabel('Amplitude (arb. units)')
    sp = sp+1;
end
legend('0.96 K','1.52 K','5.71 K','Location','northeast')


%% Plot AoT vs. temp from my MATLAB fft
for jj = 1:length(obj.FFT.range)
%     jj
    for ii = 1:length(obj.mass.range(1).upPeak)
        figure
%         subplot(2,2,ii)
        plot(obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoT,'*',...
            obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoTfitV)
%         maxFreq = obj.mass.range.upPeak(kk).maxFreq;
%         title(num2str(mean(maxFreq)))
        xlabel('Temperature (K)')
        ylabel('Amplitude/T (arb. units)')
        
    end
%     suptitle('MATLAB ft, Fourier Window 15 to 35 T')
end

%% Plot mass vs center field and frequency vs center field
% figure
for ii = 1:length(obj.FFT.range)
    CF(ii) = mean(obj.FFT.range(ii).upTemp(1).range);
    for jj = [1 2]%:length(obj.mass.range(ii).upPeak)
%         chjj = jj
        massCF(jj,ii) = obj.mass.range(ii).upPeak(jj).m;
        errorCF(jj,ii) = obj.mass.range(ii).upPeak(jj).AoTrms;
        maxPeakSTD(jj,ii) = std(obj.mass.range(ii).upPeak(jj).maxFreq);
        maxPeakCF(jj,ii) = mean(obj.mass.range(ii).upPeak(jj).maxFreq);
    end
end

figure(5)
figure(6)
for kk = [1 2] %:length(maxPeakCF(:,1))
    figure(5)
    
    errorbar(CF,maxPeakCF(kk,:),maxPeakSTD(kk,:),'*')
    hold on
    
    figure(6)
    title('center field vs mass')
    plot(CF,massCF(kk,:),'*')
%     errorbar(CF,massCF(kk,:),errorCF(kk,:),'-*')
%     legend(
    hold on
end
figure(5)
title('center field vs max peaks')
legend('Peek 1','Peek 2','Peek 4','Peek 5')
xlabel('Center Field (T)')
ylabel('Peek Frequency (T)')

figure(6)
title('center field vs mass')
legend('Peek 1','Peek 2','Peek 4','Peek 5')
xlabel('Center Field (T)')
ylabel('Effective Mass')


%% Plot max peak amplitude vs center peaks
peaks = [1 2 4 5];
figure(1)
figure(2)
for kk = 1:4
    for ii = 1:length(obj.mass.range)
        CF(ii) = mean(obj.FFT.range(ii).upTemp(1).range);
        A(:,ii) = obj.mass.range(ii).upPeak(kk).A;
        phase(:,ii) = obj.mass.range(ii).upPeak(kk).phase;
%         if isempty(Aval) == 0
%             A(:,ii) = Aval;
%         end
    end
   
    
    for jj = 1:length(A(:,1))
        figure(1)
        subplot(2,2,kk)
        leg{jj} = strcat(num2str(obj.FFT.range(kk).upTemp(jj).temp), 'K');
        plot(CF,A(jj,:))
        xlabel('Center Field (T)')
        ylabel('Fourier amp. (arb. units)')
        hold on
        
        figure(2)
        subplot(2,2,kk)
        leg{jj} = strcat(num2str(obj.FFT.range(kk).upTemp(jj).temp), 'K');
        plot(CF,phase(jj,:))
        xlabel('Center Field (T)')
        ylabel('phase (rad)')
        hold on
    end

    title(strcat('Peak',num2str(peaks(kk))))
%     leg
    for kk = 4
        figure(1)
        subplot(2,2,kk)
        legend(leg)
        
        figure(2)
        subplot(2,2,kk)
        legend(leg)
    end
%     pause(1)
end