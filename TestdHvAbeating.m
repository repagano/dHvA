close all
clear
clc
%% view interference on 
H = 10:.00001:55;
xi = cos(2*pi*3757./H+pi*3757/2)+cos(2*pi*3695./H);%
objFab.xUp = H;
objFab.yUp = xi;
objFab.temp = 1;
%%
endFields = 15:.5:55;%[35];%
icenterFields = 1./endFields;
iFFspan = abs((1/40-1/55));%.02;%
fRange = [4 4500];
for jj = 1:length(endFields)
    iendField = 1/endFields(jj);
    iFFrange = [iendField+iFFspan, iendField];
    FFrange = 1./iFFrange     ;                  
    FFT = FFTcalc(objFab,'dHvA',FFrange,'up',fRange,'filter');      
    FFTmaxFab(jj) = max(FFT.FFT);
    FFTrangeFab(jj) = mean(FFT.range);
end
% objFab.FFT = FFTcalc(objFab,endFields,iFFspan,[5 4500]);%R
FFTfab.max = FFTmaxFab;
FFTfab.peak = FFTrangeFab;
%%
figure
plot(1./objFab.xUp,objFab.yUp)
%%
figure
plot(1./FFTrangeFab,FFTmaxFab)

%% Load data and Plot max of hig pass filtered data
% close all
% clear
clc
%% Use class to load data 
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\031919_PdGa\031919\';
temps = [1.52 .62 5.71 3.96 8 10 12 14 18 3.27 2.24 .96];
files = dir(strcat(loc,'*.ASC'));
filesFFT = dir(strcat(loc,'FT*'));
[temps,I] = find(temps == .96);
files = files(I);
filesFFT = filesFFT(I);
%%
obj = dHvA(loc,temps,files,filesFFT);
%%
endFields = 15:.5:55;%
icenterFields = 1./endFields;
iFFspan = abs((1/40-1/55));
% %% Check span
% icenterField = 1./endFields;
% for ii = 1:length(endFields) 
%     chii = ii
%     s = [icenterFields(ii) icenterField(ii)+iFFspan]
%     span = 1./s
% end

%% Find maxs and plot
dHvA.FFTload(obj,endFields,iFFspan,[5 4500]);
%%
for ii = 1:length(obj.FFT.range)
    [FFTmax(ii), I] = max(obj.FFT.range(ii).upTemp.FFT);
    CFmax(ii) = mean(obj.FFT.range(ii).upTemp.range);
end

%%
figure
plot(1./CFmax,FFTmax)
hold on
% figure
plot(1./FFTfab.peak,FFTfab.max*.025/.73)
title('Beating between 3757 and 3698')
legend('FFT','simulated FFT')
xlabel('1/CF (1/T)')
ylabel('Fourier Amp. (Arb. Units)')
%%
figure
for ii = 1:2:length(obj.FFT.range)

plot(obj.FFT.range(ii).upTemp.f,obj.FFT.range(ii).upTemp.FFT)
hold on
pause(1)
end