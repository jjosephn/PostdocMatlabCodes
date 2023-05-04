function [morphFeatures] = morphologicalFeatures(FHRUCSignal,segLenParam,sampFreq,plotFlag)

if(nargin<4)
    plotFlag = 0;
end

FHRSignal = FHRUCSignal(1,:);
UCSignal = FHRUCSignal(2,:);

segLenMinitue = segLenParam.segLenMinitue;
baselineSubSegLenMinitue = segLenParam.baselineSubSegLenMinitue;
baselineSubSegSamples = baselineSubSegLenMinitue*60*sampFreq;
varSegLenMinitue = segLenParam.varSegLenMinitue;
variabilitySegmentSamples = varSegLenMinitue*60*sampFreq;
[morphFeatures.FHRBaselineSegment, continuousSegments] = FHRBaselineEstimation(FHRSignal);

[UCSignalInfo] = UCInfoExtract(UCSignal,sampFreq);

morphFeatures.UCSignalLocations = UCSignalInfo.UCSignalLocations;
morphFeatures.UCSignalPeaks = UCSignalInfo.UCSignalPeaks;
morphFeatures.UCSignalPeaksNumbers = UCSignalInfo.UCSignalPeaksNumbers;
morphFeatures.UCSignalStartLocation = UCSignalInfo.UCSignalStartLocation;
morphFeatures.UCSignalStartVal = UCSignalInfo.UCSignalStartVal;
morphFeatures.UCSignalStopLocation = UCSignalInfo.UCSignalStopLocation;
morphFeatures.UCSignalStopVal = UCSignalInfo.UCSignalStopVal;

if length(continuousSegments) >= baselineSubSegSamples
    FHRAccelerationSegmentIndex = FHRSignal >= (morphFeatures.FHRBaselineSegment + 15);

    accelerationCount = 0;
    accelerationBelow15Count = 0;
    accelerationAbove15Count = 0;
    accelFlag1 = 0;
    accelFlag2 = 0;
    for accelVar1=1:1:length(FHRAccelerationSegmentIndex)
        if FHRAccelerationSegmentIndex(accelVar1) == 1 && accelFlag1 == 0
            accelerationStartIndex = accelVar1;
            accelFlag1 = 1;
        elseif FHRAccelerationSegmentIndex(accelVar1) == 0 && accelFlag1 == 1
            accelerationStopIndex = accelVar1;
            accelFlag1 = 0;
            accelFlag2 = 1;
        end
        while accelFlag2
            accelerationCount = accelerationCount + 1;
            accelerationDuration = (accelerationStopIndex - accelerationStartIndex)/sampFreq;
            if accelerationDuration <= 15
                accelerationBelow15Count = accelerationBelow15Count + 1;
            elseif accelerationDuration > 15
                accelerationAbove15Count = accelerationAbove15Count + 1;
            end
            accelFlag2 = 0;
        end

    end

    morphFeatures.FHRAccelerationSegmentNum = accelerationCount;
    morphFeatures.FHRAccelerationSegmentDurationBelow15 = accelerationBelow15Count;
    morphFeatures.FHRAccelerationSegmentDurationAbove15 = accelerationAbove15Count;
    
    
    FHRBaselineSegment = morphFeatures.FHRBaselineSegment;
    baselineFHR = (FHRBaselineSegment)/abs(max(FHRSignal));
    baselineAbove15 = (FHRBaselineSegment + 15)/abs(max(FHRSignal));
    baselineBelow15 = (FHRBaselineSegment - 15)/abs(max(FHRSignal));
    FHRSignalPlot = FHRSignal/abs(max(FHRSignal));
    UCSignalPlot = UCSignal/abs(max(UCSignal));


    if plotFlag == 1
        figure
        plot(2+baselineAbove15*ones(1,length(FHRSignal)),'b');
        hold on;
        plot(2+baselineFHR*ones(1,length(FHRSignal)),'r');
        hold on;
        plot(2+baselineBelow15*ones(1,length(FHRSignal)),'b');
        hold on;
        plot(2+FHRSignalPlot,'b');
        hold on;
        plot(1+UCSignalPlot,'g');
        hold on;
        plot(UCSignalInfo.UCSignalLocations,1+UCSignalInfo.UCSignalPeaks,'r*');
        hold on;
        plot(UCSignalInfo.UCSignalStartLocation,1 + UCSignalInfo.UCSignalStartVal,'bo');
        hold on;
        plot(UCSignalInfo.UCSignalStopLocation,1 + UCSignalInfo.UCSignalStopVal,'ro');
        hold on;
        title('CTG Feature Extraction')
        xlabel('Samples')
        ylabel('Normalized CTG (FHR and UC)')
    end
    
    
    
    FHRDecelerationBaseLineSegmentIndex = FHRSignal <= (morphFeatures.FHRBaselineSegment);    
    decelerationBaselineStartIndex = zeros(1,length(FHRDecelerationBaseLineSegmentIndex));
    decelerationBaselineStopIndex = zeros(1,length(FHRDecelerationBaseLineSegmentIndex));

    decelStartFlag = 0;
    decelStopFlag = 1;
    decelVar2 = 1;
    for decelVar1=1:1:length(FHRDecelerationBaseLineSegmentIndex)

        if FHRDecelerationBaseLineSegmentIndex(decelVar1) == 1 && decelStartFlag == 0 && decelStopFlag == 1
            decelerationBaselineStartIndex(decelVar2) = decelVar1;
            decelStartFlag = 1;
            decelStopFlag = 0;
        end
        if FHRDecelerationBaseLineSegmentIndex(decelVar1) == 0 && decelStartFlag == 1 && decelStopFlag == 0
            decelerationBaselineStopIndex(decelVar2) = decelVar1;
            decelVar2 = decelVar2 + 1;
            decelStartFlag = 0;
            decelStopFlag = 1;
        elseif decelStopFlag == 0 && decelStartFlag == 1 && decelVar1 == length(FHRDecelerationBaseLineSegmentIndex)
            decelerationBaselineStopIndex(decelVar2) = decelVar1;
            decelVar2 = decelVar2 + 1;
            decelStartFlag = 0;
            decelStopFlag = 1;
        end
    end
    
    decelerationBaselineStartIndex(decelerationBaselineStartIndex==0) = [];
    decelerationBaselineStopIndex(decelerationBaselineStopIndex==0) = [];
        
    decelerationStartIndex = zeros(1,length(decelerationBaselineStartIndex));
    decelerationStopIndex = zeros(1,length(decelerationBaselineStopIndex));

    
    for decelVar1=1:1:length(decelerationBaselineStartIndex)
        currentDecelerationSegmentIndex = FHRSignal(decelerationBaselineStartIndex(decelVar1):decelerationBaselineStopIndex(decelVar1))  <= (morphFeatures.FHRBaselineSegment - 15);
        if length(find(currentDecelerationSegmentIndex>0))/sampFreq >= 15
            decelerationStartIndex(decelVar1) = decelerationBaselineStartIndex(decelVar1);
            decelerationStopIndex(decelVar1) = decelerationBaselineStopIndex(decelVar1);
        end
    end
   
    decelerationStartIndex(decelerationStartIndex==0) = [];
    decelerationStopIndex(decelerationStopIndex==0) = [];
    
    
    if plotFlag == 1
        plot(decelerationStartIndex,2 + FHRSignalPlot(decelerationStartIndex),'g*');
        hold on;
        plot(decelerationStopIndex,2 + FHRSignalPlot(decelerationStopIndex),'r*');
        hold on
    end
    
    decelerationCount = 0;
    decelerationDuration = 0;
    earlyDecelerationCount =  0;
    variableDecelerationCount =  0;
    lateDecelerationCount = 0;
    prolongedDecelerationCount = 0;

    
    for decelVar1=1:1:length(decelerationStartIndex)
        [~, nadirLocation] = min(FHRSignal(decelerationStartIndex(decelVar1):decelerationStopIndex(decelVar1)));
        nadirLocation = decelerationStartIndex(decelVar1) + nadirLocation - 1;
        startToNadir = (nadirLocation - decelerationStartIndex(decelVar1))/sampFreq;
        nadirToStop = (decelerationStopIndex(decelVar1) - nadirLocation)/sampFreq;
        
        decelerationCount = decelerationCount + 1;
        decelerationDuration = decelerationDuration + (decelerationStopIndex(decelVar1) - decelerationStartIndex(decelVar1))/sampFreq;
        
        if startToNadir < 30
            if startToNadir <= 15
                earlyDecelerationCount = earlyDecelerationCount + 1;
            else
                variableDecelerationCount = variableDecelerationCount + 1;
            end
        end
        
         if startToNadir >= 30 || nadirToStop >=30
            lateDecelerationCount = lateDecelerationCount + 1;
         end
        
         if ((decelerationStopIndex(decelVar1) - decelerationStartIndex(decelVar1))/sampFreq) >= 180
            prolongedDecelerationCount = prolongedDecelerationCount + 1;
        end
    end
    
    morphFeatures.FHRDecelerationSegmentNum = decelerationCount;
    morphFeatures.FHRDecelerationSegmentDuration = decelerationDuration;
    morphFeatures.FHRDecelerationSegmentEarlyDecelerationNum = earlyDecelerationCount;
    morphFeatures.FHRDecelerationSegmentVariableDecelerationNum = variableDecelerationCount;
    morphFeatures.FHRDecelerationSegmentLateDecelerationNum = lateDecelerationCount;
    morphFeatures.FHRDecelerationSegmentProlongedDecelerationNum = prolongedDecelerationCount;


    morphFeatures.FHRSegmentVariability = zeros(1,segLenMinitue/varSegLenMinitue);
    startVar = length(FHRSignal) - variabilitySegmentSamples + 1;
    stopVar = length(FHRSignal(1,:));
    flag5 = 1;
    var4 = 1;
    while flag5
        morphFeatures.FHRSegmentVariability(var4) = max(FHRSignal(1,startVar:stopVar)) - min(FHRSignal(1,startVar:stopVar));
        stopVar = startVar - 1;
        startVar = stopVar - variabilitySegmentSamples + 1;
        var4 = var4 + 1;
        if startVar <= 0 && stopVar > 0
            morphFeatures.FHRSegmentVariability(var4) = max(FHRSignal(1,1:stopVar)) - min(FHRSignal(1,1:stopVar));
            flag5 = 0;
        elseif stopVar <= 0
            flag5 = 0;
        end
    end
else
    morphFeatures.FHRAccelerationSegmentNum = 0;
    morphFeatures.FHRAccelerationSegmentDurationBelow15 = 0;
    morphFeatures.FHRAccelerationSegmentDurationAbove15 = 0;
    morphFeatures.FHRDecelerationSegmentNum = 0;
    morphFeatures.FHRDecelerationSegmentDuration = 0;
    morphFeatures.FHRDecelerationSegmentEarlyDecelerationNum = 0;
    morphFeatures.FHRDecelerationSegmentVariableDecelerationNum = 0;
    morphFeatures.FHRDecelerationSegmentLateDecelerationNum = 0;
    morphFeatures.FHRDecelerationSegmentProlongedDecelerationNum = 0;
    morphFeatures.FHRSegmentVariability = zeros(1,segLenMinitue/varSegLenMinitue);
end
