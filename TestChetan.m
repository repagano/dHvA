close all
clear 
clc
%%
filename = 'C:\Users\repag\Documents\DiTusa\MagLab\NALGERESISTANCE.037\';%NALGERESISTANCE.025.txt';
files = dir(strcat(filename,'*.txt'));
files = files(16:end);
for ii = 1:length(files)
    delimiter = '\t';
    startRow = 8;
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%*s%[^\n\r]';
    filesFFT = strcat(filename,files(ii).name);
    names{ii} = files(ii).name;
    fileID = fopen(filesFFT,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false);
    fclose(fileID);

    H = dataArray{1}; %field in Tesla    
    v = dataArray{6}; %voltage for resistance
    c = dataArray{8}; %capacitance in pico farad
    raw(ii).H = H;
    raw(ii).v = v;
    raw(ii).c = c;
    
    if (H(1)-H(end))>0
        H = flip(H);
        v = flip(v);
        c = flip(c);
    end

    objFabV.xUp = H;
    objFabV.yUp = v;
    objFabV.temp = 1;
    FFrange = [H(30),H(end)];
    fRange = [1 100];
%     FFTv(ii) = FFTcalc(objFabV,'whatever',FFrange,'up',fRange,'filter'); 
    
    
    objFabC.xUp = H;
    objFabC.yUp = c;
    objFabC.temp = 1;
    FFTc(ii) = FFTcalc(objFabC,'whatever',FFrange,'up',fRange,'filter'); 
    
    objFabCd.xUp = H(1:end-1);
    objFabCd.yUp = diff(c);
    objFabCd.temp = 1;
%     FFTcd(ii) = FFTcalc(objFabCd,'whatever',FFrange,'up',fRange,'filter');
   
    
end


% %%
% figure
% for ii = 1:length(FFTc)
%     
%     plot(FFTv(ii).f,FFTv(ii).FFT)
%     hold on
% end
% title('FFT of voltage')
% xlabel('1/H')
% ylabel('v')

%%
ii = 0
Ft1(ii,:)=FFTc(1).FFT;
ft1(ii,:) = FFTc(1).f;
Ft11(ii,:)=FFTc(11).FFT;
ft11(ii,:)=FFTc(11).f;
Ft12(ii,:)=FFTc(12).FFT;
ft12(ii,:)=FFTc(12).f;
Ft13(ii,:)=FFTc(13).FFT;
ft13(ii,:)=FFTc(13).f;
Ft14(ii,:)=FFTc(14).FFT;
ft14(ii,:)=FFTc(14).f;
Ft15(ii,:)=FFTc(15).FFT;
ft15(ii,:)=FFTc(15).f;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:11


polyV = [0:10];
[~,I] = max(F1(ii,:));
FFTmax1(ii) = f1(I);
[~,I] = max(F11(ii,:));
FFTmax11(ii) = f11(I);
[~,I] = max(F12(ii,:));
FFTmax12(ii) = f12(I);
[~,I] = max(F13(ii,:));
FFTmax13(ii) = f13(I);
[~,I] = max(F14(ii,:));
FFTmax14(ii) = f14(I);
[~,I] = max(F15(ii,:));
FFTmax15(ii) = f15(I);


figure(2)
plot(f1(ii,:),F1(ii,:))
title(names(10))
hold on

figure(3)
plot(FFTc(11).f,F11(ii,:))
title(names(11))
hold on

figure(4)
plot(FFTc(12).f,F12(ii,:))
title(names(12))
hold on

figure(5)
plot(FFTc(13).f,F13(ii,:))
title(names(13))
hold on

figure(6)
plot(FFTc(14).f,F14(ii,:))
title(names(14))
hold on

figure(7)
plot(FFTc(15).f,F15(ii,:))
title(names(15))
hold on
end
%%
leg = {'poly 0 fit','poly 1 fit','poly 2 fit','poly 3 fit','poly 4 fit',...
    'poly 5 fit','poly 6 fit','poly 7 fit','poly 8 fit','poly 9 fit','poly 10 fit'}
figure(2)
legend(leg)

figure(3)
legend(leg)

figure(4)
legend(leg)

figure(5)
legend(leg)

figure(6)
legend(leg)

figure(7)
legend(leg)
%%
figure(1)
plot(polyV,FFTmax1,'*-',polyV,FFTmax12,'*-',polyV,FFTmax13,'*-',...
    polyV,FFTmax14,'*-',polyV,FFTmax15,'*-')
hold on
legend(names)
xlabel('degree of polynomial fit')
ylabel('frequency of peak')
%%
% figure(1)
% figure(2)
% for ii = 1:length(FFTcd)
% %     figure
% %     plot(FFTcd(ii).f,FFTcd(ii).FFT)%
% %     title(strcat('line',num2str(ii))) 
% %     hold on
%     figure
%     plot(FFTc(ii).f,FFTc(ii).FFT)
%     title(strcat('line',num2str(ii))) 
% %     hold on
% %     pause(1)
% end
% figure(1)
% title('FFT of dC vs 1/H degree 10')
% xlabel('frequency')
% ylabel('amplitude')
% 
% figure(2)
% title('FFT of c vs 1/H polinomial degree 5')
% xlabel('frequency')
% ylabel('amplitude')
%%
% 
for ii = 11;%1:length(FFTc)
    figure
    plot(FFTc(ii).x,FFTc(ii).y,'r')%,FFTc(ii).xspl,FFTc(ii).yspl,'.b')
    title(strcat(names(ii),'no fit'))
    xlabel('1/H')
    ylabel('c')
%     hold on
end
%%
figure
plot(FFTc(3).x,FFTc(3).y)

figure
plot(FFTc(11).x,FFTc(11).y)
%%
figure
plot(raw(15).H,raw(15).c)

figure
plot(raw(14).H,raw(14).c)

figure
plot(raw(11).H,raw(11).c)

figure
plot(raw(3).H,raw(3).c)
% 
% for ii = 1:length(FFTcd)
%     figure
%     plot(FFTcd(ii).x,FFTcd(ii).y,'r',FFTcd(ii).xspl,FFTcd(ii).yspl,'.b')
%     hold on
% end
% %%
% 
% xinv = objFabV.xUp;
% FF = 1./xinv;
% xi = objFabV.yUp;
% FFInd = FF >= FFrange(1) & FF <= FFrange(2);
% x = FF(FFInd);
% y = xi(FFInd);
% 
% figure
% plot(x,y)
% 
% n = 2^15;
% xspl = linspace(x(end)