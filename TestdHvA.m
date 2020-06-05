clear
close all
clc
 %% Here the files to be loaded are chosen.
addpath C:\Users\repag\Documents\MATLAB\DiTusa
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\031919_PdGa\031919\';
temps = [1.52 .62 5.71 3.96 8 10 12 14 18 3.27 2.24 .96];%Temperatures
% corresponding to files to be loaded
files = dir(strcat(loc,'*.ASC'));
filesFFT = dir(strcat(loc,'FT*'));
% here set the variable varTemps equal to the temperatures you wish to
% load. If you set varTemps = temps all files are loaded.
% [~,I] = max(temps)
varTemps = temps;
I = ismember(temps,varTemps);
temps = temps(I);
files = files(I);
filesFFT = filesFFT(I);

%% Here the raw data is loaded into dHvA class
obj = dHvA(loc,temps,files,filesFFT);

%% Here the 1/H windows that are selected and fourier transformed are defined
endFields = 35;%the endFields are the maximum field values of each window
iFFspan = abs((1/15-1/35));%Here the width of the 1/H window is defined 

%% Here the Fourier transform is performed  
dHvA.FFTload(obj,endFields,iFFspan,[5 4500],'extFFT');

%% Here the effective masses of the different peaks are calculated 
% peakRange is an matrix whose rows define the range around each peak 
% identified on the Fourier spectrum. The columns correspond to different
% peaks. The range should be large enough to allow for variation in the position 
% of the peak as the different ranges and temperatures are scanned over, but 
% not be so large that another peak is selected. 
peakRange = [490 708; 708 850;  2600 2900; 3546 3880]; 
% dHvA.massLoad runs the method that calculates the effective masses.
dHvA.massLoad(obj,peakRange);

%% Plot phase for 4 different peaks, all temps 

figure
for ii = 1:4
    obj.mass.range.upPeak(ii).phase([2 9 10 11 12])...
        = obj.mass.range.upPeak(ii).phase([2 9 10 11 12])+pi
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).phase)
    hold on
%     pause(1)
end
legend('peak 1','peak 1','peak 3','peak 4')
xlabel('Temp')
ylabel('phase')
%% Phase and amplitude subplot
% for ii = 1:length
Fmax = [obj.mass.range.upPeak(1).maxFreq(1),obj.mass.range.upPeak(2).maxFreq(1),....
    obj.mass.range.upPeak(3).maxFreq(1),obj.mass.range.upPeak(4).maxFreq(1)];

Pmax = [obj.mass.range.upPeak(1).phase(1),obj.mass.range.upPeak(2).phase(1),....
    obj.mass.range.upPeak(3).phase(1),obj.mass.range.upPeak(4).phase(1)];

figure
subplot(2,1,1)
plot(obj.FFT.range.upTemp(1).f,obj.FFT.range.upTemp(1).phase,Fmax,Pmax,'*r')
title('phase from FFT')
subplot(2,1,2)
plot(obj.FFT.range.upTemp(1).f,obj.FFT.range.upTemp(1).FFT)
title('amplitude from FFT')

% figure
% plot(obj.FFT.range.upTemp(1).xspl,obj.FFT.range.upTemp(1).yspl)
%% Plot phase for 4 different peaks, all temps 
% 
% figure
% for ii = 1:4
%     ii
%     p=obj.mass.range.upPeak(ii).phase
%     I = p<pi
%     p(I) = p(I)+pi
%     phase(:,ii) = p
%     plot(obj.mass.range.upPeak(ii).temp,phase(:,ii))
%     hold on
% end
%% plot dM/dH vs 1/H. Fig 1 in document created using this section with
% endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
temp = 1:11;%[1];
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
%used to make Fig.1
%% Fourier transform from LANL and MATLAB. Fig 2 in document created using 
% this section with endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
tempFFT = 1:11;
pos = 1;

for ii = 1;%:12
%     subplot(2,2,pos)
    figure
    plot(obj.raw(ii).f,obj.raw(ii).FFT,'LineWidth',1.5)
    hold on
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).FFT)
    pos = pos+1;
%     for jj = 1:4
%         plot(obj.mass.raw(jj).maxFreq(ii),obj.mass.raw(jj).A(ii),'*')
%     end
end
title('Fig 2: 0.96 K, fourier window 15 to 35 T')
legend('LANL fourier transform','MATLAB fourier transform')
xlabel('dHvA frequency, F (T)')
ylabel('Amplitude (arb. units)')
%used to make Fig.2.1
%%
figure
for ii = 1;%:12
    subplot(1,2,1)
    plot(obj.raw(ii).f,obj.raw(ii).FFT,'LineWidth',1.5)
    hold on
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).FFT)
    pos = pos+1;
%     for jj = 1:4
%         plot(obj.mass.raw(jj).maxFreq(ii),obj.mass.raw(jj).A(ii),'*')
%     end
end
title('Fig 2.1U: 0.96 K, fourier window 15 to 35 T')
legend('LANL fourier transform','MATLAB fourier transform')
xlabel('dHvA frequency, F (T)')
ylabel('Amplitude (arb. units)')

subplot(1,2,2)
for ii =1:4
    plot(obj.mass.raw(ii).temp,obj.mass.raw(ii).A,'LineWidth',2)
    hold on
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).A)
    
end
title('Fig 2.2U')
xlabel('Temp (K)')
ylabel('Fourier Amplitude (arb. units)')
legend('Frequency ~600 matlab','Frequency ~600 Lanl','Frequency ~775 matlab','Frequency ~775 Lanl',...
    'Frequency ~2770 matlab','Frequency ~2770 Lanl','Frequency ~3700 matlab','Frequency ~3700 Lanl')
%% Compleat phase plot
% this section with endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
% tempFFT = 1:11;
pos = 1;
    
for ii = 1:length(obj.FFT.range(1).upTemp)
%     obj.FFT.range.upTemp(ii).phase(1:2:end) = obj.FFT.range.upTemp(ii).phase(1:2:end)+pi
%     subplot(2,2,pos)
    figure
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).phase)
    pos = pos+1;
   pause(1)

%     for jj = 1:4
%         plot(obj.mass.raw(jj).maxFreq(ii),obj.mass.raw(jj).A(ii),'*')
%     end
end
% title('Fig 2: 0.96 K, fourier window 15 to 35 T')
% legend('LANL fourier transform','MATLAB fourier transform')
% xlabel('dHvA frequency, F (T)')
% ylabel('Amplitude (arb. units)')


%% Plot AoT vs. temp from my MATLAB fft
for jj = 1:length(obj.FFT.range)
%     jj
figure
    for ii = 1:length(obj.mass.range(1).upPeak)
%         figure(jj)
        subplot(2,2,ii)
        plot(obj.mass.range(jj).upPeak(ii).temp,(obj.mass.range(jj).upPeak(ii).AoT),'*',...
            obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoTfitV)
%         maxFreq = obj.mass.range.upPeak(kk).maxFreq;
%         title(num2str(mean(maxFreq)))
        maxFreq = obj.mass.range(jj).upPeak(ii).maxFreq;
        title(strcat('Peak Frequency (T) =',num2str(round(mean(maxFreq)))));
        xlabel('Temperature (K)')
        ylabel('Amplitude/T (arb. units)')
        
    end
    suptitle('MATLAB ft, Fourier Window 15 to 35 T')
end
%used to make fig. 3.a
%% AoT vs temp plot from LANL fft
figure
for kk = 1:4
    subplot(2,2,kk)
    plot(obj.mass.raw(kk).temp,obj.mass.raw(kk).AoT,'*',...
        obj.mass.raw(kk).temp,obj.mass.raw(kk).AoTfitV)
    maxFreq = obj.mass.raw(kk).maxFreq;
    title(strcat('Peak Frequency (T) =',num2str(round(mean(maxFreq)))));
    xlabel('Temperature (K)')
    ylabel('Amplitude/T (arb. units)')
    
end

suptitle('LANL ft, Fouriet Window 15 to 35 T')
%used to make fig. 3.b
%% Plot masses vs frequency
figure

for jj = 1:length(obj.FFT.range)
%     jj
    for ii = 1:length(obj.mass.range(1).upPeak)
        massPlot(ii) = obj.mass.range(jj).upPeak(ii).m;
        freqPlot(ii) = mean(obj.mass.range(jj).upPeak(ii).maxFreq);
        errorPlot(ii) = mean(obj.mass.range(jj).upPeak(ii).AoTrms);
        rawMassPlot(ii) = obj.mass.raw(ii).m;
        rawFreqPlot(ii) = mean(obj.mass.raw(ii).maxFreq);
        rawErrorPlot(ii) = mean(obj.mass.raw(ii).AoTrms);
    end
    leg{jj} = num2str(obj.FFT.range(jj).upTemp(1).range);
%     errorbar(freqPlot,massPlot,errorPlot,'o')
    plot(freqPlot,massPlot,'o',rawFreqPlot,rawMassPlot,'*')
    hold on

end
% errorbar(rawFreqPlot,rawMassPlot,rawErrorPlot,'*')
title('Fourier Window 15 to 35')
xlabel('Effective mass, m*')
ylabel('F(T)')
legend('m* from MATLAB ft','m* from LANL ft')
%used to make fig.4
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
    ii
    
    subplot(2,3,sp)    
    plot(obj.FFT.range(ii).upTemp(2).f,obj.FFT.range(ii).upTemp(2).FFT,'r')
    hold on
    plot(obj.FFT.range(ii).upTemp(3).f,obj.FFT.range(ii).upTemp(3).FFT,'b')
    plot(obj.FFT.range(ii).upTemp(7).f,obj.FFT.range(ii).upTemp(7).FFT,'Color',[0 .7 .3])
    titl1 = strcat('CF = ',num2str(round(mean(obj.FFT.range(ii).upTemp(1).range))))
    titl = strcat(titl1,'T')
    title(titl)
    xlabel('frequency (T)')
    ylabel('Amplitude (arb. units)')
    sp = sp+1
end
legend('0.96 K','1.52 K','5.71 K','Location','northeast')




    

%% Plot mass vs center field and frequency vs center field
% figure
for ii = 1:length(obj.FFT.range)
    CF(ii) = mean(obj.FFT.range(ii).upTemp(1).range);
    for jj = [1:4]%:length(obj.mass.range(ii).upPeak)
%         chjj = jj
        massCF(jj,ii) = obj.mass.range(ii).upPeak(jj).m;
        errorCF(jj,ii) = obj.mass.range(ii).upPeak(jj).AoTrms;
        maxPeakSTD(jj,ii) = std(obj.mass.range(ii).upPeak(jj).maxFreq);
        maxPeakCF(jj,ii) = mean(obj.mass.range(ii).upPeak(jj).maxFreq);
    end
end

figure(5)
figure(6)
for kk = [1 4] %:length(maxPeakCF(:,1))
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
%% calcluate and plot masses multipl times to test fit
% figure
% for jj = 1:15
%     chjj = jj
%     dHvA.massLoad(obj,peakRange)
%     for ii = 1:length(obj.mass.range(1).upPeak)
%         plot(mean(obj.mass.range(1).upPeak(ii).maxFreq),obj.mass.range(1).upPeak(ii).m,'*')%,'o',obj.mass.raw(ii).maxFreq,obj.mass.raw(ii).m,'*')%,obj.mass.down(ii).maxFreq,obj.mass.down(ii).m,'x',
%         hold on
%     end
% pause(.5)
% end
% xlabel('Frequency (T)')
% ylabel('mass')
% legend('my FT','LANL FT')


