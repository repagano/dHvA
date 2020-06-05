clear
close all
clc
%% Here the files to be loaded are chosen.
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\031919_PdGa\031919\';
temps = [1.52 .62 5.71 3.96 8 10 12 14 18 3.27 2.24 .96];%Temperatures
% corresponding to files to be loaded
files = dir(strcat(loc,'*.ASC'));
filesFFT = dir(strcat(loc,'FT*'));
% here set the variable temps equal to the data temperatures corresponding 
% to the files you wish to load in the expression find(temps == temps) by 
% setting temps == temps, all files are loaded.
[~,I] = find(temps == temps);
temps = temps(I);
files = files(I);
filesFFT = filesFFT(I);

%% Here the raw data is loaded into dHvA class
obj = dHvA(loc,temps,files,filesFFT);

%% Here the 1/H windows that are selected and fourier transformed are defined
endFields = [15:10:55];%the endFields are the maximum field values of each window
iFFspan = abs((1/40-1/55));%Here the width of the 1/H window is defined 

%% Here the Fourier transform is performed  
dHvA.FFTload(obj,endFields,iFFspan,[5 4500]);

%% Here the effective masses of the different peaks are calculated 
% peakRange is an matrix whose rows define the range around each peak 
% identified on the Fourier spectrum. The columns correspond to different
% peaks. The range should be large enough to allow for variation in the position 
% of the peak as the different ranges and temperatures are scanned over, but 
% not be so large that another peak is selected. 
peakRange = [490 708; 708 850;  2600 2900; 3546 3880]; 
% dHvA.massLoad runs the method that calculates the effective masses.
dHvA.massLoad(obj,peakRange);

%% plot dM/dH vs 1/H. Fig 1 in document created using this section with
% endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
temp = 1:length(obj.raw);
for jj = length(obj.FFT.range)
    for ii = temp
        
        figure(ii)
        subplot(2,1,1)
        plot(obj.FFT.range(jj).upTemp(ii).x,obj.FFT.range(jj).upTemp(ii).y,'Linewidth',1.2)
        hold on
        plot(obj.FFT.range(jj).upTemp(ii).xspl,obj.FFT.range(jj).upTemp(ii).yspl,'Linewidth',.4)
        xlabel('1/H')
        ylabel('dM/dH')
        legend('dM/dH vs 1/H','Spline fit')
        subplot(2,1,2)
        plot(1./obj.FFT.range(jj).upTemp(ii).x,obj.FFT.range(jj).upTemp(ii).y,'Linewidth',1)
        hold on
        plot(1./obj.FFT.range(jj).upTemp(ii).xspl,obj.FFT.range(jj).upTemp(ii).yspl,'Linewidth',.4)
        xlabel('H')
        ylabel('dM/dH')
        legend('dM/dH vs 1/H','Spline fit')
    end      
end
suptitl1 = strcat(num2str(obj.raw(ii).temp),'K');
suptitle(strcat('Fig 1: T=',suptitl1))
%% plot Fourier transform for all temps 
temp = [2];% 3 4 5 6 7 10];
sp = 1;
figure
leg = [];
for jj = 1:length(obj.FFT.range)
%     subplot(2,3,sp)
    for ii = temp%1:length(obj.raw)
%         figure(jj)
%         leg{jj} = strcat(CFnum2str(round(mean(obj.FFT.range(jj).upTemp(ii).range))),'T');
obj.FFT.range(jj).upTemp(ii).range
        txt = {sprintf('CF = %.0f T',round(mean(obj.FFT.range(jj).upTemp(ii).range)))};
        leg = [leg,txt];
%         subplot(2,4,sp)
        plot(obj.FFT.range(jj).upTemp(ii).f,obj.FFT.range(jj).upTemp(ii).FFT)
        hold on
        
    end

sp = sp+1;
% titl1 = strcat('CF = ',num2str(round(mean(obj.FFT.range(jj).upTemp(1).range))));
% titl = strcat(titl1,'T');
txt = sprintf('1/H window = %0.1d (1/T)',iFFspan)
title(txt)  
% title(titl)
xlabel('frequency (T)')
ylabel('Amplitude (arb. units)')

end
legend(leg)
% Used to creat fig.5a-c (fig 5c has a different lable but is the third figure) 
% allthought two different ranges was scanned over. 
%% Plot max peaks vs. temps
clearvars leg
for jj = 1:length(obj.mass.range)
    %%
    figure
%     jj = 8;
    for ii = 1:length(obj.mass.range(jj).upPeak)
        plot(obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).A,'-*')
        hold on
        leg{ii} = strcat('f = ', num2str(round(mean(obj.mass.range(jj).upPeak(ii).maxFreq))));
    end
    xlabel('T (K)')
    ylabel('Amplitude (arb. units)')
    title(num2str(obj.FFT.range(jj).upTemp(1).range))
    legend(leg)
end
%% Plot AoT vs. temp from my MATLAB fft

for jj = 1:length(obj.FFT.range)
    %%
    figure
%     jj = 8;
    for ii = 1:length(obj.mass.range(1).upPeak)
        
        subplot(2,2,ii)
        plot(obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoT,'*',...
            obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoTfitV)
%         maxFreq = obj.mass.range.upPeak(kk).maxFreq;
%         title(num2str(mean(maxFreq)))
        xlabel('Temperature (K)')
        ylabel('Amplitude/T (arb. units)')
        
    end
    titl = num2str(mean(obj.FFT.range(jj).upTemp(1).range));
    suptitle(titl)
end
endCnt = 23
%% Plot mass vs center field and frequency vs center field
% figure
if ~exist('endCnt')
    endCnt = 1;
end
% chendCnt = endCnt
for ii = 1:length(obj.FFT.range)
    CF(ii) = mean(obj.FFT.range(ii).upTemp(1).range);
    for jj = [1:4]%:length(obj.mass.range(ii).upPeak)
%         chjj = jj
        massCF(jj,ii) = obj.mass.range(ii).upPeak(jj).m;
        errorCF(jj,ii) = obj.mass.range(ii).upPeak(jj).AoTrms;
        maxPeakSTD(jj,ii) = std(obj.mass.range(ii).upPeak(jj).maxFreq);
        maxPeakCF(jj,ii) = mean(obj.mass.range(ii).upPeak(jj).maxFreq);
        maxPhaseCF(jj,ii) = mean(obj.mass.range(ii).upPeak(jj).phase);
    end
end

figure(endCnt)
figure(endCnt+1)
figure(endCnt+2)

for kk = 1:length(maxPeakCF(:,1))
    figure(endCnt)
    errorbar(CF,maxPeakCF(kk,:),maxPeakSTD(kk,:),'-*')
    hold on
    
    figure(endCnt+1)
    plot(CF,massCF(kk,:),'-*')
%     errorbar(CF,massCF(kk,:),errorCF(kk,:),'-*')
%     legend(
    hold on
    
    figure(endCnt+2)
%     p = maxPhaseCF<pi;
%     maxPhaseCF(I) = maxPhaseCF(I)+pi;
    plot(CF,maxPhaseCF(kk,:),'*-')
    hold on
%     pp = pp+1
end
figure(endCnt)
title('center field vs max peaks')
legend('Peak 1','Peak 2','Peak 4','Peak 5')
xlabel('Center Field (T)')
ylabel('Peak Frequency (T)')

figure(endCnt+1)
title('Fig. 6.b m* vs center field')
legend('Peak 1','Peak 2','Peak 4','Peak 5')
xlabel('Center Field (T)')
ylabel('effective mass')

figure(endCnt+2)
title('phase vs center field')
legend('Peak 1','Peak 2','Peak 4','Peak 5')
xlabel('Center field (T)')
ylabel('phase')


% Used to creat Fig 6 allthought two different ranges was scanned over.