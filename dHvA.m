classdef dHvA < handle
    properties
        dataType        
        raw
        FFT
        mass
        scan
    end

    methods(Static)
%%
        function obj=dHvA(loc,temps,varargin)            
        % dHvA initializes the class and loads the data to be analyzed.
        %
        % obj = dHvA(loc, temps, files, filesFFT) initializes the class obj
        % and defines the property raw. loc is the string indicating the 
        % file location, temps is the array of temperatures corresponding  
        % to each file, files is the struct array of files containing the 
        % raw measured data, and filesFFT is the struct array the files 
        % containing previously calculated FFT data. dHvA defines the 
        % property raw: a struct array containing  xUp which is the rising 
        % external field, yUp which is the magnetic susceptibility of the 
        % material as the external field rises; xDown and yDown are
        % respectively the external field and magnetic susceptibility as 
        % the field falls; temp is the temperature at which the data was 
        % taken; f and FFT is the 
        % frequency and Fourier amplitude of the previously calculated FFT; 
        % range is the range of the previously calculated FFT assumed to be 
        % 15 to 35 T. The property dataType indicates what type of data was 
        % loaded: de Haas van Alphan (dhva) measurement or a rf penetration
        % depth (rfPD) measurement (NOTE, the program is only designed to 
        % load dhva data when a previously calculated FFT is present. The 
        % dhva data is identified by a .ASC file).
        %
        % obj = dHvA(loc, temp, files) is used when there is not previously
        % calculated FFT.  The function will
        % assume dhva data if the file type is .ASC and rfPD data if the 
        % file type is .dat. No other file types are supported. 

            fclose all;
            [sortedTemp,sortTempI] = sort(temps);
            names = varargin{1};
            
            for ii = 1:length(names) 
                % Organize files in order of ascending temperature and 
                % create string of file location and name 
                I = sortTempI(ii);
                name = names(I).name;
                file = strcat(loc,name);
                
                
                if contains(name,'.ASC') && length(varargin) == 2
                    % Load dhva data and FFT data
                    namesFFT = varargin{2};
                    nameFFT = namesFFT(I).name;
                    fileFFT = strcat(loc,nameFFT);
                    raw(ii) = dHvALoad(sortedTemp(ii),file,fileFFT);
                    dataType = 'dHvA';
                elseif contains(name,'.ASC') && length(varargin) == 1
                    % Load dhva data
                    raw(ii) = dHvALoad(sortedTemp(ii),file);
                    dataType = 'dHvA';
                elseif contains(name,'.dat')
                    % Load rfPD data
                    raw(ii) = PenetrationDepthLoad(sortedTemp(ii),file);
                    dataType = 'rfPD';
                else
                    error('Inappropriate inputs to function')
                end
            end
            % Define properties raw and dataType
            obj.raw = raw;
            obj.dataType = dataType;
        end
       
%%        
        function FFTload(obj,endFields,iFFspan,varargin)
        % FFTload selects the input H range, Fourier transforms that range,
        % then selects input frequency range.
        
        % dHvA.FFTload(obj,endFields,iFFspan,[fRangeTop fRangeBottom]) will 
        % select the data falling inside of the specified H range window and
        % Fourier transformed it. The window is define by the 1/H window
        % (iFFspan) and the strength of the field that defines the end of
        % the range (endFields). If an array is input instead of a single
        % value here (endFields), multiple ranges will be scanned over. The
        % data is then input to the function FFTcalc and Fourier transformed. 
        
        % dHvA.FFTload(obj,endFields,iFFspan,[fRangeTop fRangeBottom], 'extFFT')
        % will Fourier transform data as described above, then rescale the
        % previously calculated FFT to match the arbitrary unites produced
        % by the function FFTcalc. 

            
            fRange = varargin{1}; 
            for jj = 1:length(endFields)
                % Calculate H range to be fourier transformed 
                iendField = 1/endFields(jj);
                iFFrange = [iendField+iFFspan iendField];
                FFrange = 1./iFFrange;
                %Fourier transform range
                for ii = 1:length(obj.raw)                           
                    upTemp(ii) = FFTcalc(obj.raw(ii),obj.dataType,...
                        FFrange,'up',fRange);
                end                 
                obj.FFT.range(jj).upTemp = upTemp;
            end
            
%           % recale previously calculated Fourier transform
            if length(varargin) == 2
                for kk = 1:length(obj.raw)
                    fraw = obj.raw(kk).f;
                    frawI = fraw>=fRange(1) & fraw<=fRange(2);
                    obj.raw(kk).f = fraw(frawI);
                    obj.raw(kk).FFT = obj.raw(kk).FFT(frawI).*9.67e6;%.*8.2e7;
                end
            end  
        end 
%%
        function massLoad(obj,peakRange)
            % massLoad calculates the effective masses of each Fourier 
            % amplitude peak identified by the user.
            %
            % dHvA.massLoad(obj, peakRange) scans through each peak 
            % identified by the user, and calculates the effective mass
            % using the function massCalc. The peaks (peakRange) are indicated using a
            % N by 2 array where N is the number of peaks. The first and
            % second columb identify the begining and end of the
            % range of interest in the data corisponding to each peak. If
            % dataType is type 'dHvA' then the masses for the function will
            % run massCalc for externally calculated data as well. 
            
            
            % massCalc is run for each peak of internally calculated peak
            for jj = 1:length(obj.FFT.range)
                cnt = 1;
                for ii = 1:length(peakRange(:,1))
                    massUpPeak(ii) = massCalc(obj.FFT.range(jj).upTemp,peakRange(cnt,:),'on');
                    cnt = cnt+1;
                end
            obj.mass.range(jj).upPeak = massUpPeak;
            end

                
            % massCalc is run for each peak of externally calculated peak
            if strcmp(obj.dataType,'dHvA')
                
                cnt = 1;
                for ii = 1:length(peakRange(:,1))
                    massRaw(ii) = massCalc(obj.raw,peakRange(cnt,:),'off');
                    cnt = cnt+1;

                end
                obj.mass.raw = massRaw;
            end
        end
    end   
end
