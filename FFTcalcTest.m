function [FFT] = FFTcalcTest(obj,dataType,FFrange,FFdir,varargin)
   % FFT = FFTcalc(obj,
%     obj
    fRange = varargin{1};
    FFdirNum = strcmp(FFdir,'up');
    Tbot = FFrange(1);
    Ttop = FFrange(2);      
    switch FFdir
        case 'up'  
            FF = obj.xUp;
            xi = obj.yUp;
%             FF(high:end) = zeros(length(FF(high:end)),1);  
            FFInd = FF >= Tbot & FF <= Ttop;
            xinv = FF(FFInd);
            int = .0000001;%change to the power of 8
            x = xinv;%1./
            N = 2^20;%round(abs(x(1)-x(end))/int); %Total number of points of array fed into FFT 
            
            if rem(N,2) == 1
                N = N+1;
            end
            
        case 'down'
            FF = obj.xDown;
            xi = obj.yDown;
            FF(1:high) = zeros(high,1);
            FFInd = FF >= Tbot & FF <= Ttop;
            xinv = FF(FFInd);
            int = .0000001;
            x = 1./xinv;
            if rem(N,2) == 1
                N = N+1;
            end
%             chN = N
            xspl = linspace(x(end),x(1),N);
    end
%     chN = N;
    n = 10000; %round(N*2/3); %number of points fit to the spline %;
    xspl = linspace(x(1),x(end),n);
    int = xspl(2)-xspl(1);
    y = xi(FFInd);
    yog = y;
    if strcmp(dataType,'dHvA')
        degree = 3;
    elseif strcmp(dataType,'rfPD')
        degree = 2;
    elseif strcmp(dataType,'whatever')
        degree = 1;
    end
    yfit = polyfit(xinv,y,degree);
    yfitV = polyval(yfit,xinv);
    y = y-yfitV;
    
    yspl = spline(x,y,xspl(1:n));
    Datapoints = length(yspl);
    hwind = hann(Datapoints)';
    
    nZero = (N-n)/2;
    xfft1 = xspl(1)-floor(nZero)*int;
    xfftEnd = xspl(end)+ceil(nZero)*int;  
    

    zeros1 = zeros(1,floor(nZero));
    zerosEnd = zeros(1,ceil(nZero));
    
    ysplW = yspl.*hwind;
    xfft = linspace(xfft1,xfftEnd,N);%xspl;%
    yfft = [zeros1,ysplW,zerosEnd];%ysplW;%
    
%     figure
%     plot(xfft,yfft)
    
    chxfft = length(xfft);
    chyfft = length(yfft);

    Fs = (abs(xfft(1)-xfft(end))/N)^-1;

    
    on = 0;
    if on == 1
        figure
        subplot(2,1,1)
        plot(xinv,yog,'b',xinv,yfitV,'r')
         subplot(2,1,2)
        plot(xinv,y)
        temp = strcat(num2str(obj.temp),'K');
        suptitle(strcat('T=', temp))
    end
    
    FFT.x = x;
    FFT.y = y;
    FFT.xspl = xspl;
    FFT.yspl = yspl;
    Datapoints = length(yfft);
    Length=abs(xfft(end)-xfft(1));

    FFTc = fft(yfft);%./length(yfft);
    FFTphase = angle(FFTc);
    FFTphaseVal = FFTphase(1:Datapoints/2+1);
    FFTval = abs(FFTc(1:Datapoints/2+1));

    fs=Datapoints/Length;  
    f = fs/2*linspace(0,1,Datapoints/2+1);
    fI = f>=fRange(1) & f<=fRange(2);
    FFT.f = f(fI);
    FFT.FFT = FFTval(fI);
    FFT.phase = FFTphaseVal(fI);
    FFT.temp = obj.temp;
    FFT.range = FFrange;

end
    


%     chFS = (abs(xspl(1)-xspl(end))/n)^-1
%     format short
    %     if length(varargin) == 2
%         chHP = 1
%         yspl = highpass(yspl,3500,Fs);
%     end


%     chN = N;
%     format short
%     y = 
    
%     figure
%     plot(xfft,yfft)
%     hold on
%     plot(xspl,ysplW)
    
%     format long



%     xZeros1 = linspace(xfft1,xspl(1)-int,floor(nZero));
%     chz1 = xZeros1(2)-xZeros1(1)
%     ch1 = length(xZeros1(1,:))
%     xZerosEnd = linspace(xfftEnd,xspl(end)+int,ceil(nZero));
%     chzend = xZerosEnd(2)-xZerosEnd(1);
%     ch2 = length(xZerosEnd(1,:))