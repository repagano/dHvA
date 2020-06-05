function raw = dHvALoad(sortedTemp, varargin)
    % raw = dHvALoad(sortedTemp, file, fileFFT) where sortedTemp is the 
    % temp associated with the file, file is the directory of the raw data 
    % file, and fileFFT is the directory of the Fourier transform data, 
    % returns struct raw with variables xUp which is the rising external 
    % field, yUp which is the magnetic susceptibility of the material as 
    % the external field rises; xDown and yDown are respectively the 
    % external field and magnetic susceptibility as the field falls; temp 
    % is the temperature at which the data was taken; f and FFT is the 
    % frequency and Fourier amplitude of the previously calculated FFT; 
    % range is the range of the previously calculated FFT assumed to be 15 
    % to 35 T.
    %
    % raw = dHvALoad(sortedTemp, file) returns the struct raw with the 
    % variables described above, excluding the files f, FFT, and range. 
    
    % Load dHvA raw data file 
    filename = varargin{1};
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
        'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines',...
        startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    % seperate data into arrays  
    Bdot = (dataArray{1});%dH/dt 
    Mag1 = (dataArray{2});%dM/dt
    FF = (dataArray{4});%fixed field
   
    % Separate data into rising and falling field data and calculate 
    % magnetic susceptibility xi
    [~,Ttop] = max(FF);
    Tbot = 1.5;
    xi = Mag1./Bdot;                
    FFup = [FF(1:Ttop);zeros(length(FF(Ttop:end)),1)];
    
    FFupInd = FFup >= Tbot & FFup <= Ttop;
    raw.xUp = FFup(FFupInd);
    raw.yUp = xi(FFupInd);

    FFdown = [zeros(Ttop,1);FF(Ttop:end)];
    FFdownInd = FFdown >= Tbot & FFdown <= Ttop;
    raw.xDown = FF(FFdownInd);
    raw.yDown = xi(FFdownInd);                
    raw.temp = sortedTemp;
    raw.dataType = 'dHvA';
    
    % Load FFT data if fileFFT exists 
    if length(varargin) == 2
        FFTfilename = varargin{2};
        FFTraw = importdata(FFTfilename);
        raw.f = FFTraw(:,1);
        raw.FFT = FFTraw(:,2);
        raw.range = [15 35];
    end
end