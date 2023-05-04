clc;
close all;
clearvars; 

[parameters,signal,info] = initDatabase();

fileVals = info.fileVals;
fileValsChar = info.fileValsChar;
parameterTypeNames = info.parameterTypeNames;
parameterNames = info.parameterNames;

pH690count = 0;
pH700To709count = 0;
pH710To719count = 0;
pH720count = 0;

pHLowCount = 0;
pHHighCount = 0;

pHLowfileVals = zeros(1,20);
pHHighfileVals = zeros(1,20);

pH690fileVals = zeros(1,20);
pH700To709fileVals = zeros(1,32);
pH710To719fileVals = zeros(1,108);
pH720fileVals = zeros(1,375);

missingSamplesFHR = zeros(1,length(info.fileVals));
missingSamplesUC = zeros(1,length(info.fileVals));

dataPercent = 0.8;
fileNumbLow = 194;    % pH <=7.2
fileNumbHigh = fileNumbLow;


countLowTrain = 0;
countLowTest = 0;
countHighTrain = 0;
countHighTest = 0;

trainFlag = 0;
testFlag = 0;
trainCount = 1;
testCount = 1;
featureFlag = 0;

for fileNo=1:1:length(info.fileVals)
    
    fprintf('Processing File No. %d - %s \n',fileNo,fileValsChar{fileNo});
    
    signalLength = parameters.(fileValsChar{fileNo}).(parameterTypeNames{1}).samplesPerSignal;
    secondStageStart = parameters.(fileValsChar{fileNo}).(parameterTypeNames{7}).PosIISt;    
    timeScale = signal.(fileValsChar{fileNo}).timeScale(1:end);
    
    pH = parameters.(fileValsChar{fileNo}).(parameterTypeNames{2}).pH;
    BDecf = parameters.(fileValsChar{fileNo}).(parameterTypeNames{2}).BDecf;
    BE = parameters.(fileValsChar{fileNo}).(parameterTypeNames{2}).BE;
    sampFreq = parameters.(fileValsChar{fileNo}).(parameterTypeNames{1}).samplingFreq;
    
    
    FHR = signal.(fileValsChar{fileNo}).FHR(1:end);    
    UC = signal.(fileValsChar{fileNo}).UC(1:end);
    
    if pH <= 7.2
        pHLowCount = pHLowCount + 1;
        pHLowfileVals(pHLowCount) = fileVals(fileNo);
        if secondStageStart > -1
            FHR1 = signal.(fileValsChar{fileNo}).FHR(1:secondStageStart-1); % Stage !
            UC1 = signal.(fileValsChar{fileNo}).UC(1:secondStageStart-1);
            FHR2 = signal.(fileValsChar{fileNo}).FHR(secondStageStart:end); % Stage 2
            UC2 = signal.(fileValsChar{fileNo}).UC(secondStageStart:end);
        else
            FHR1 = signal.(fileValsChar{fileNo}).FHR(1:end);
            UC1 = signal.(fileValsChar{fileNo}).UC(1:end);
        end    
        if pHLowCount <= floor(dataPercent*fileNumbLow)
            trainFlag = 1;
            featureFlag = 1;
            countLowTrain = countLowTrain + 1;
        elseif pHLowCount > floor(dataPercent*fileNumbLow) && pHLowCount <= fileNumbLow
            testFlag = 1;
            featureFlag = 1;
            countLowTest = countLowTest + 1;
        end
    elseif pH > 7.2
        pHHighCount = pHHighCount + 1;
        pHHighfileVals(pHHighCount) = fileVals(fileNo);
        if secondStageStart > -1
            FHR1 = signal.(fileValsChar{fileNo}).FHR(1:secondStageStart-1); % Stage 1
            UC1 = signal.(fileValsChar{fileNo}).UC(1:secondStageStart-1);
            FHR2 = signal.(fileValsChar{fileNo}).FHR(secondStageStart:end); % Stage 2
            UC2 = signal.(fileValsChar{fileNo}).UC(secondStageStart:end);
        else
            FHR1 = signal.(fileValsChar{fileNo}).FHR(1:end);
            UC1 = signal.(fileValsChar{fileNo}).UC(1:end);
        end        
        if pHHighCount <= floor(dataPercent*fileNumbHigh)
            trainFlag = 1;
            featureFlag = 1;
            countHighTrain = countHighTrain + 1;
        elseif pHHighCount > floor(dataPercent*fileNumbHigh) && pHHighCount <= fileNumbHigh
            testFlag = 1;
            featureFlag = 1;
            countHighTest = countHighTest + 1;
        end
    end
    
    
    
    %%

    while ((trainFlag || testFlag) && featureFlag)
        [missingSamplesFHRStart,missingSamplesFHRStop] = missingSamples(FHR);  % Location of missing samples
        [missingSamplesUCStart,missingSamplesUCStop] = missingSamples(UC);

        segLenSeconds = 20;             % in seconds according to FIGO
        skipSegmentSamples = segLenSeconds*sampFreq;

        tempFullSignalStage1 = [FHR1;UC1];
        if secondStageStart > -1
            tempFullSignalStage2 = [FHR2;UC2];
        else
            tempFullSignalStage2 = zeros(2,1);
        end
        selectFlag = 1; % 0 removes missing samples more than 20 seconds in FHR only. 1 removes missing samples more than 20 seconds from both FHR and UC.
        rmMissing = removeMissingSamples(tempFullSignalStage1,tempFullSignalStage2,selectFlag,sampFreq,secondStageStart); %remove missing samples more than 20 seconds
        tempFullSignalStage11 = rmMissing.tempFullSignalStage11;  % 11 Stage1 FHR
        tempFullSignalStage12 = rmMissing.tempFullSignalStage12;    % 12 Stage 1 UC
        tempFullSignalStage21 = rmMissing.tempFullSignalStage21;    % Stage 2 FHR
        tempFullSignalStage22 = rmMissing.tempFullSignalStage22;    % Stage 2 UC

        interMissing = interpolateMissingSamples(rmMissing,secondStageStart);

        fullSignalStageI = interMissing.fullSignalStageI;
        fullSignalStageII = interMissing.fullSignalStageII;
        
        if isempty(fullSignalStageI)
            trainFlag = 0;
            testFlag = 0;
        else

            segLenMinitue = 10;                 % Change to thedesired segment length in minutes. Should not be less than 10. FIGO guidelines for ref.
            baselineSegmentSamples = segLenMinitue*60*sampFreq;
            slopeSegLenMinitue = .5;
            slopeSegSamples = slopeSegLenMinitue*60*sampFreq;
            baselineSubSegLenMinitue = 2;
            baselineSubSegSamples = baselineSubSegLenMinitue*60*sampFreq;
            varSegLenMinitue = 1;
            variabilitySegmentSamples = varSegLenMinitue*60*sampFreq;

            segLenParam.segLenMinitue = segLenMinitue;
            segLenParam.slopeSegLenMinitue = slopeSegLenMinitue;
            segLenParam.baselineSubSegLenMinitue = baselineSubSegLenMinitue;
            segLenParam.varSegLenMinitue = varSegLenMinitue;

            statFeatMean = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatMedian = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatSD = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatSE = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatVar = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatMax = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatMin = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatMaxMinDiff = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatMeanAbsDev = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   
            statFeatMedianAbsDev = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));   

            baselineSegmentFHR = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));    
            accelerationNoSegmentFHR = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            accelerationDurationBelow15SegmentFHR = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            accelerationDurationAbove15SegmentFHR = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            decelerationNoSegmentFHR = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            decelerationDuration = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            earlyDecelerationCount = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            variableDecelerationCount = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            lateDecelerationCount = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            prolongedDecelerationCount = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            meanVariabilitySegmentFHR = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            varVariabilitySegmentFHR = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));
            variabilitySegmentFHR = zeros(segLenMinitue/varSegLenMinitue,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));

            UCSignalPeaksNumbers = zeros(1,ceil(length(fullSignalStageI(1,:))/baselineSegmentSamples));

            start = length(fullSignalStageI(1,:)) - baselineSegmentSamples + 1;
            stop = length(fullSignalStageI(1,:));
            flag1 = 1;
            var1 = 1;
            while flag1
                if start <= 1
                    signalFHRUC = fullSignalStageI(:,1:stop);
                    flag1 = 0;
                else
                    signalFHRUC = fullSignalStageI(:,start:stop);
                end

                statFeatures = statisticalFeatures(signalFHRUC(1,:));

                statFeatMean(var1) = statFeatures.statFeatMean;
                statFeatMedian(var1) = statFeatures.statFeatMedian;
                statFeatSD(var1) = statFeatures.statFeatSD;
                statFeatSE(var1) = statFeatures.statFeatSE;
                statFeatVar(var1) = statFeatures.statFeatVar;
                statFeatMax(var1) = statFeatures.statFeatMax;
                statFeatMin(var1) = statFeatures.statFeatMin;
                statFeatMaxMinDiff(var1) = statFeatures.statFeatMaxMinDiff;
                statFeatMeanAbsDev(var1) = statFeatures.statFeatMeanAbsDev;
                statFeatMedianAbsDev(var1) = statFeatures.statFeatMedianAbsDev;

                morphFeatures = morphologicalFeatures(signalFHRUC,segLenParam,sampFreq);

                baselineSegmentFHR(var1) = morphFeatures.FHRBaselineSegment;

                accelerationNoSegmentFHR(var1) = morphFeatures.FHRAccelerationSegmentNum;
                accelerationDurationBelow15SegmentFHR(var1) = morphFeatures.FHRAccelerationSegmentDurationBelow15;
                accelerationDurationAbove15SegmentFHR(var1) = morphFeatures.FHRAccelerationSegmentDurationAbove15;

                decelerationNoSegmentFHR(var1) = morphFeatures.FHRDecelerationSegmentNum;
                decelerationDuration(var1) = morphFeatures.FHRDecelerationSegmentDuration;
                earlyDecelerationCount(var1) = morphFeatures.FHRDecelerationSegmentEarlyDecelerationNum;
                variableDecelerationCount(var1) = morphFeatures.FHRDecelerationSegmentVariableDecelerationNum;
                lateDecelerationCount(var1) = morphFeatures.FHRDecelerationSegmentLateDecelerationNum;
                prolongedDecelerationCount(var1) = morphFeatures.FHRDecelerationSegmentProlongedDecelerationNum;
                variabilitySegmentFHR(:,var1) = morphFeatures.FHRSegmentVariability';
                meanVariabilitySegmentFHR(var1) = mean(variabilitySegmentFHR(:,var1));
                varVariabilitySegmentFHR(var1) = var(variabilitySegmentFHR(:,var1));

                UCSignalPeaksNumbers(var1) = morphFeatures.UCSignalPeaksNumbers;

                stop = start - 1;
                start = stop - baselineSegmentSamples + 1;
                var1 = var1 + 1;        
            end
            featureFlag = 0;
        end
    end
    
    while trainFlag
        features.train.fileInfo{trainCount} = fileValsChar{fileNo};
        features.train.(fileValsChar{fileNo}).clinicalParamsVectorTrain = [pH BDecf BE];
        features.train.(fileValsChar{fileNo}).featureVectorTrain = [statFeatMean; statFeatMedian; statFeatSD; statFeatSE; statFeatVar; statFeatMax; ... 
            statFeatMin; statFeatMaxMinDiff; statFeatMeanAbsDev; statFeatMedianAbsDev; baselineSegmentFHR; ... 
            accelerationNoSegmentFHR; accelerationDurationBelow15SegmentFHR; accelerationDurationAbove15SegmentFHR; ... 
            decelerationNoSegmentFHR; decelerationDuration; earlyDecelerationCount; variableDecelerationCount; ... 
            lateDecelerationCount; prolongedDecelerationCount; meanVariabilitySegmentFHR; varVariabilitySegmentFHR;
            UCSignalPeaksNumbers]; 
%             variabilitySegmentFHR];
        trainCount = trainCount + 1;
        trainFlag = 0;
    end
    
    while testFlag
        features.test.fileInfo{testCount} = fileValsChar{fileNo};
        features.test.(fileValsChar{fileNo}).clinicalParamsVectorTest = [pH BDecf BE];  
        features.test.(fileValsChar{fileNo}).featureVectorTest = [statFeatMean; statFeatMedian; statFeatSD; statFeatSE; statFeatVar; statFeatMax; ... 
            statFeatMin; statFeatMaxMinDiff; statFeatMeanAbsDev; statFeatMedianAbsDev; baselineSegmentFHR; ... 
            accelerationNoSegmentFHR; accelerationDurationBelow15SegmentFHR; accelerationDurationAbove15SegmentFHR; ... 
            decelerationNoSegmentFHR; decelerationDuration; earlyDecelerationCount; variableDecelerationCount; ... 
            lateDecelerationCount; prolongedDecelerationCount; meanVariabilitySegmentFHR; varVariabilitySegmentFHR;
            UCSignalPeaksNumbers];
%             variabilitySegmentFHR];
        testCount = testCount + 1;  
        testFlag = 0;        
    end 
    clear fullSignalStageIDetrend fullSignalStageIIDetrend;
end

% savefile = 'featuresExtractionV1o010MtsSegmentEqual';
% % savefile = 'featuresExtractionV1o010MtsSegment';
% save(savefile,'features');
