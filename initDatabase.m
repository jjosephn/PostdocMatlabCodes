function [parameters,signal,info] = initDatabase()

% Read data from the database. parameters contains data from .hea file
% (clinical and signal parameters). signal contains the fetal heart rate
% (FHR) and uterine contractions (UC). Info contains file numbers, parameter names etc.

addpath('... /CTGDatabase/...');                                            % add location to the database


fileVals = [1001:1506 2001:2046];
fileValsChar{length(fileVals)} = [];
parameterTypeNames{7} = [];
parameterTypeNames{1} = 'SignalParameters';
parametersTypeCount = 1;

parameterNames{35} = [];

for fileNo = 1:1:length(fileVals)
        fileValsChar{fileNo} = strcat('file',int2str(fileVals(fileNo)));
        info.headerFileName{fileNo} = strcat(int2str(fileVals(fileNo)),'.hea');
        info.dataFileName{fileNo} = strcat(int2str(fileVals(fileNo)),'.dat');


        fid1 = fopen(info.headerFileName{fileNo});
        tline = fgetl(fid1);
        tlineData = sscanf(tline,'%d %d %d %d');
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).recordName = num2str(tlineData(1));
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).noOfSignals = tlineData(2);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).samplingFreq = tlineData(3);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).samplesPerSignal = tlineData(4);

        clear tlineData

        tline = fgetl(fid1);
        tlineData = sscanf(tline,'%*s %d %d(0)/bpm %d %d %d %d %d %*s');
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).signalFormatFHR = tlineData(1);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCGainFHR = tlineData(2);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCResolutionFHR = tlineData(3);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCOffsetFHR = tlineData(4);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).initialValueFHR = tlineData(5);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).checksumFHR = tlineData(6);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).blockSizeFHR = tlineData(7);

        clear tlineData

        tline = fgetl(fid1);
        tlineData = sscanf(tline,'%*s %d %d/nd %d %d %d %d %d %*s');
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).signalFormatUC = tlineData(1);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCGainUC = tlineData(2);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCResolutionUC = tlineData(3);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCOffsetUC = tlineData(4);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).initialValueUC = tlineData(5);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).checksumUC = tlineData(6);
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).blockSizeUC = tlineData(7);

        parametersNameCount = 1;
        eof = 0;
        while ~eof
            tline = fgetl(fid1);
            if tline == -1
                eof = 1;
            elseif ~isempty(tline) && strcmp(tline(1:4),'#-- ')                
                parametersTypeCount = parametersTypeCount + 1;
                parameterTypeName = strtrim(tline(5:end));
                parameterTypeNames{parametersTypeCount} = renameFunction(parameterTypeName);
            elseif ~isempty(tline) && ~strcmp(tline(1),'#-- ') && ~strcmp(tline(1:4),'#---')
                parameterName = strtrim(tline(2:14));
                parameterNames{parametersNameCount} = renameFunction(parameterName);
                parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).(parameterNames{parametersNameCount}) = str2double(tline(15:end));     
                parametersNameCount = parametersNameCount + 1;
            end
        end

        fclose(fid1);


        % read data from a file
        fid1 =fopen(info.dataFileName{fileNo},'r');
        parametersTypeCount = 1;
        rawData = fread(fid1, [parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).noOfSignals, ... 
            parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).samplesPerSignal],'uint16', 0);
        fclose(fid1);

        signal.(fileValsChar{fileNo}).FHR = (rawData(1,:) - parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCOffsetFHR)/... 
        parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCGainFHR;
        signal.(fileValsChar{fileNo}).UC = (rawData(2,:) - parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCOffsetUC)/...
            parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).ADCGainUC;
        signal.(fileValsChar{fileNo}).timeScale = (0:(parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).samplesPerSignal - 1))/...
            parameters.(fileValsChar{fileNo}).(parameterTypeNames{parametersTypeCount}).samplingFreq;


end

info.fileValsChar = fileValsChar;
info.fileVals = fileVals;
info.parameterTypeNames = parameterTypeNames;
info.parameterNames = parameterNames;
info.fileVals = fileVals;



