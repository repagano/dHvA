close all
clc
%%
f = obj.FFT.range.upTemp(2).f';
FFT = obj.FFT.range.upTemp(2).FFT';
data1 = [f,FFT];
save('FT35to55_temp0.96','data1','-ascii')

f = obj.FFT.range.upTemp(3).f';
FFT = obj.FFT.range.upTemp(3).FFT';
data1 = [f,FFT];
save('FT35to55_temp1.52','data1','-ascii')

f = obj.FFT.range.upTemp(7).f';
FFT = obj.FFT.range.upTemp(7).FFT';
data1 = [f,FFT];
save('FT35to55_temp5.71','data1','-ascii')

%%

L1 = load('FT35to55_temp0.96')%'FT25.67to35_temp0.96')%
L2 = load('FT35to55_temp1.52')%'FT25.67to35_temp1.52')%
L3 = load('FT35to55_temp5.71')%'FT25.67to35_temp5.71')%
% L1 = ans;

figure
plot(L1(:,1),L1(:,2),L2(:,1),L2(:,2),L3(:,1),L3(:,2))

%%

L1 = load('FT25.67to35_temp0.96')%
L2 = load('FT25.67to35_temp1.52')%
L3 = load('FT25.67to35_temp5.71')%
% L1 = ans;

figure
plot(L1(:,1),L1(:,2),L2(:,1),L2(:,2),L3(:,1),L3(:,2))