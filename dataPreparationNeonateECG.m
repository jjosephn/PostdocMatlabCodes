clc;
close all;
clearvars; 
addpath('.../NeonateECGDataCork/mat_files');                                % Add path to the neonate data

cellPatientInfo = table2cell(readtable('patientInfo.xls'));
cellParam = table2cell(readtable('timeVaryingParam.xls'));

[totalRecords,~] = size(cellPatientInfo);


segmentNum = 1;
for recordNo = 1:1:totalRecords
    
    filename = strcat('ID',string(cellPatientInfo{recordNo,1}));
    currentFileName = strcat(filename,'.mat');
    recordName = sprintf('ID%dpart%d',cellPatientInfo{recordNo,1},0);
    if exist(currentFileName, 'file') == 2
        currentNeonateData.(recordName) = load(currentFileName);
    end

    var1 = 1;
    while var1 <= 7
        recordName = sprintf('ID%dpart%d',cellPatientInfo{recordNo,1},var1);
        currentFileName = strcat(filename,sprintf('_part%d.mat',var1));
        if exist(currentFileName, 'file') == 2
            currentNeonateData.(recordName) = load(currentFileName);
        end
        var1 = var1 + 1;
    end

    recordNameAll = fieldnames(currentNeonateData);
    for recordNum = 1:1:length(fieldnames(currentNeonateData))
        record = currentNeonateData.(recordNameAll{recordNum}).data;
        sampFreq = currentNeonateData.(recordNameAll{recordNum}).Fs;
        channelLabels = currentNeonateData.(recordNameAll{recordNum}).ch_labels;
        [noOfChannesls, noOfSamples] = size(record);
        ecgStartPNA = currentNeonateData.(recordNameAll{recordNum}).start_ECG_pna_hours;    
        currentRecordNum = cellPatientInfo{recordNo,1};

        preProcessedSignal = record;
        [~,sampPoints] = size(cellParam);
        refNum = find([cellParam{:,1}] == currentRecordNum);
        currentRecordParam = cellParam(refNum,2:sampPoints);    

        ecgStopPNA = ecgStartPNA + noOfSamples/(sampFreq*60*60);                % in hours
        segmentMinutes = 10;                                                    % Change to the required time. 10 minutes in this case.
        segmentLength = segmentMinutes*60*sampFreq;
        a = 1;
        for var1=2:1:length(currentRecordParam)
            sampleTimeMinutes = currentRecordParam{1,var1};
            sampleTimeHours = sampleTimeMinutes/60;

            if sampleTimeHours >= ecgStartPNA && sampleTimeHours <= ecgStopPNA
                segmentName = sprintf('segment%d',segmentNum);
                samplePoint = ceil((sampleTimeHours - ecgStartPNA)*3600*sampFreq);
                segmentStart = samplePoint - segmentLength + 1;
                segmentStop = samplePoint + segmentLength;
                if segmentStart <= 0
                    segmentStart = 1;
                end            
                if segmentStop > length(preProcessedSignal(1,:))
                    segmentStop = length(preProcessedSignal(1,:));
                end
                
                
                segPreProcessedSignal = zeros(2,segmentStop - segmentStart + 1);
                if strcmp(channelLabels{1},'EKG L') == 1 || strcmp(channelLabels{1},'ECG1') == 1
                    neonateECGSegment.(segmentName).channel(1) = 1;
                    segPreProcessedSignal(1,:) = preProcessedSignal(1,segmentStart:segmentStop);
                elseif strcmp(channelLabels{1},'EKG R') == 1 || strcmp(channelLabels{1},'ECG2') == 1
                    neonateECGSegment.(segmentName).channel(1) = 2;
                    segPreProcessedSignal(2,:) = preProcessedSignal(1,segmentStart:segmentStop);
                elseif strcmp(channelLabels{1},'EKG') == 1
                    neonateECGSegment.(segmentName).channel(1) = 1;
                    segPreProcessedSignal(1,:) = preProcessedSignal(1,segmentStart:segmentStop);                    
                    segPreProcessedSignal(2,:) = preProcessedSignal(1,segmentStart:segmentStop);
                end
                if length(channelLabels) > 1
                    if strcmp(channelLabels{2},'EKG L') == 1 || strcmp(channelLabels{2},'ECG1') == 1
                        neonateECGSegment.(segmentName).channel(2) = 1;
                        segPreProcessedSignal(1,:) = preProcessedSignal(2,segmentStart:segmentStop);
                    elseif strcmp(channelLabels{2},'EKG R') == 1 || strcmp(channelLabels{2},'ECG2') == 1
                        neonateECGSegment.(segmentName).channel(2) = 2;
                        segPreProcessedSignal(2,:) = preProcessedSignal(2,segmentStart:segmentStop);
                    end
                end
                
                neonateECGSegment.(segmentName).channelName = channelLabels;                
                neonateECGSegment.(segmentName).segmentSamplePoint = segmentLength;
                neonateECGSegment.(segmentName).record = segPreProcessedSignal;
                neonateECGSegment.(segmentName).sampFreq = sampFreq;
                neonateECGSegment.(segmentName).recordName = recordNameAll{recordNum};
                neonateECGSegment.(segmentName).param.pH = currentRecordParam{2,var1};
                neonateECGSegment.(segmentName).param.pCO2 = currentRecordParam{3,var1};
                neonateECGSegment.(segmentName).param.pO2 = currentRecordParam{4,var1};
                neonateECGSegment.(segmentName).param.Glucose = currentRecordParam{5,var1};
                neonateECGSegment.(segmentName).param.Lactate = currentRecordParam{6,var1};
                neonateECGSegment.(segmentName).param.Bilirubin = currentRecordParam{7,var1};
                neonateECGSegment.(segmentName).param.Base = currentRecordParam{8,var1};
                neonateECGSegment.(segmentName).param.HCO3 = currentRecordParam{9,var1};

                if length(channelLabels) > 1
                    fprintf('SegmentNo:%d; Channel Name: %s,%s; Channel No: %d,%d; \n',segmentNum,channelLabels{1}, ... 
                        channelLabels{2},neonateECGSegment.(segmentName).channel(1),neonateECGSegment.(segmentName).channel(2));
                else
                    fprintf('SegmentNo:%d; Channel Name: %s; Channel No: %d; \n',segmentNum,channelLabels{1}, ...
                       neonateECGSegment.(segmentName).channel(1));
                end
                segmentNum = segmentNum + 1;                    
                clear segPreProcessedSignal;
            end
        end
    end
    clear currentNeonateData;
end

save('neonateECGSegment600Sec','neonateECGSegment');                            % save to filename.mat



