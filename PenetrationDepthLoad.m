function raw = PenetrationDepthLoad(sortedTemp,filename)
    % raw = PenetrationDepthLoad(sortedTemp, file, fileFFT) where 
    % sortedTemp is the temp associated with the file, file is the 
    % directory of the raw data file, and fileFFT is the directory of the 
    % Fourier transform data, returns struct raw with variables xUp which 
    % is the rising external field, yUp which is the magnetic 
    % susceptibility of the material as the external field rises; xDown and 
    % yDown are respectively the external field and magnetic susceptibility 
    % as the field falls; temp is the temperature at which the data was 
    % taken.
    
    % Load rfPD raw data file 
    delimiter = '\t';
    formatSpec = '%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
        'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
    fclose(fileID);
    % seperate data into arrays  
    FF = dataArray{1}; %fixed field
    delF = dataArray{2}; 
    
    % Separate data into rising and falling field data
    [~,Ttop] = max(FF);
    Tbot = 1.5;               
    FFup = [FF(1:Ttop);zeros(length(FF(Ttop:end)),1)];
    FFupInd = FFup >= Tbot & FFup <= Ttop;
    raw.xUp = FFup(FFupInd);
    raw.yUp = delF(FFupInd);

    FFdown = [zeros(Ttop,1);FF(Ttop:end)];
    FFdownInd = FFdown >= Tbot & FFdown <= Ttop;
    raw.xDown = FF(FFdownInd);
    raw.yDown = delF(FFdownInd);                
    raw.temp = sortedTemp;
%     raw.dataType = 'rfDT';
end