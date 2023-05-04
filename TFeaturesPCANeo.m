function [TFeatures] = TFeaturesPCANeo(QRSFeatures,sampFreq)



signalBandPass = QRSFeatures.signalBandPass;
peaksLocationBandPass = QRSFeatures.peaksLocationBandPass;

inSignal = signalBandPass;
inSignal = inSignal - mean(inSignal);
inSignal = inSignal/max(abs(inSignal));

%% %%  PCA T-wave limits estimation

RPeaksIndex = peaksLocationBandPass;

RRSamples = zeros(1,length(RPeaksIndex)-1);
for var1=1:1:length(RPeaksIndex)-1
    RRSamples(var1) = RPeaksIndex(var1 + 1) - RPeaksIndex(var1);
end


maxSamples = sampFreq;
splitRecord = zeros(length(RRSamples) - 1,maxSamples);
count = 1;
for RRNum=1:1:length(RRSamples)
    beginSample = RPeaksIndex(RRNum) - ceil(RRSamples(RRNum)*1/5);
    endSample = RPeaksIndex(RRNum) + ceil(RRSamples(RRNum)*4/5);
    if (beginSample > 0)
        windowedSignal = inSignal(beginSample:endSample);
        windowedSignalSequence = 1:length(windowedSignal);
        sampleWidth = length(windowedSignal)/maxSamples;
        sampleTime = 0:sampleWidth:length(windowedSignal)-sampleWidth;
        splitRecord(count,:) = spline(windowedSignalSequence,windowedSignal,sampleTime);
        count = count + 1;
        clear windowedSignal windowedSignalSequence sampleWidth sampleTime;
    end
end

splitRecord = splitRecord';

covSplitRecord = cov(splitRecord);
[eigVecSplitRecord, eigValSplitRecord] = eig(covSplitRecord,'nobalance');
[~,orderSplitRecord] = sort(diag(eigValSplitRecord),'descend');  %# sort eigenvalues in descending order
eigVecSplitRecord = eigVecSplitRecord(:,orderSplitRecord);

neonateECGPCA = splitRecord*eigVecSplitRecord;

neonateECGPCA = neonateECGPCA';

neonateECGPromPCA = neonateECGPCA(1,:);
neonateECGPromPCA = neonateECGPromPCA - mean(neonateECGPromPCA);
neonateECGPromPCA = neonateECGPromPCA/max(abs(neonateECGPromPCA));


derivativeTWaveSignal = diff(neonateECGPromPCA);
derivativeTWaveSignal = [derivativeTWaveSignal 0];

derivativeTWaveSignal = derivativeTWaveSignal - mean(derivativeTWaveSignal);
derivativeTWaveSignal = derivativeTWaveSignal/max(abs(derivativeTWaveSignal));

[~, neonateECGPromPCAPeakLoc] = findpeaks(neonateECGPromPCA);

neonateECGPromPCAPeakLocUpdated = neonateECGPromPCAPeakLoc(neonateECGPromPCAPeakLoc<round(maxSamples*2/5));
meanVal = mean(neonateECGPromPCA(neonateECGPromPCAPeakLocUpdated));
neonateECGPromPCAPeakLocUpdated2 = neonateECGPromPCAPeakLocUpdated(neonateECGPromPCA(neonateECGPromPCAPeakLocUpdated)>=meanVal);

if length(neonateECGPromPCAPeakLocUpdated2) > 1
    [~,maxLoc] = max(neonateECGPromPCA(neonateECGPromPCAPeakLocUpdated2));
    neonateECGPromPCAPeakLocUpdated2 = neonateECGPromPCAPeakLocUpdated2(maxLoc); 
end



windoWBeginSample = round(0.2*maxSamples);
windowMiddleSample1 = round(0.30*maxSamples);
windowMiddleSample2 = round(0.50*maxSamples);
windowEndSample = round(0.65*maxSamples);


currentWindow = zeros(1,maxSamples);
windowStart = neonateECGPromPCAPeakLocUpdated2;

for windowSampleNum = windowStart + windoWBeginSample:1:windowStart + windowMiddleSample1
    currentWindow(windowSampleNum) = ((windowSampleNum  - windowStart - windoWBeginSample)/(windowMiddleSample1 - windoWBeginSample));
end
currentWindow(windowStart + windowMiddleSample1 + 1:windowStart + windowMiddleSample2) = 1;
for windowSampleNum = windowStart + windowMiddleSample2 + 1:1:windowStart + windowEndSample
    currentWindow(windowSampleNum) = 1 - ((windowSampleNum - windowStart - windowMiddleSample2)/(windowEndSample - windowMiddleSample2));
end


windowWidth = find(currentWindow>0);
winDerivativeTWaveSignal = derivativeTWaveSignal.*currentWindow;

[maxValue,maxValueLocation] = max(winDerivativeTWaveSignal(windowStart + windoWBeginSample:windowStart + windowEndSample));
[minValue,minValueLocation] = min(winDerivativeTWaveSignal(windowStart + windoWBeginSample:windowStart + windowEndSample));

maxValueLocation = windowStart + windoWBeginSample + maxValueLocation - 1;
minValueLocation = windowStart + windoWBeginSample + minValueLocation - 1;

if maxValueLocation < minValueLocation
    TWaveShapePromPCA = 1;
    thresholdTWave = minValue/2;
    
    var3 = minValueLocation;
    var4 = 1;
    thresholdFlag = 1;
    while thresholdFlag
        currentValue = winDerivativeTWaveSignal(var3);
        if var3 == windowWidth(length(windowWidth))
            TOffLocationsPromPCA = minValueLocation + var4 - 1;
            thresholdFlag = 0;
        elseif currentValue > thresholdTWave
            TOffLocationsPromPCA = minValueLocation + var4 - 1;
            thresholdFlag = 0;
        end
        var3 = var3 + 1;
        var4 = var4 + 1;
    end

    diffValue = round(0.10*(TOffLocationsPromPCA - (windowStart + windoWBeginSample)));
    [~,temp1] = max(neonateECGPromPCA(windowStart + windoWBeginSample:TOffLocationsPromPCA - diffValue));
    TPeakLocationsPromPCA = windowStart + windoWBeginSample + temp1 - 1;

    TPeakValuePromPCA = neonateECGPromPCA(TPeakLocationsPromPCA);            
    TOffValuePromPCA  = neonateECGPromPCA(TOffLocationsPromPCA);

elseif maxValueLocation > minValueLocation        
    TWaveShapePromPCA = -1;
    thresholdTWave = maxValue/2;     

    var3 = maxValueLocation;
    var4 = 1;
    thresholdFlag = 1;
    while thresholdFlag
        currentValue = winDerivativeTWaveSignal(var3);
        if var3 == windowWidth(length(windowWidth))
            TOffLocationsPromPCA = maxValueLocation + var4 - 1;
            thresholdFlag = 0;
        elseif currentValue < thresholdTWave
            TOffLocationsPromPCA = maxValueLocation + var4 - 1;
            thresholdFlag = 0;
        end
        var3 = var3 + 1;
        var4 = var4 + 1;
    end

    diffValue = round(0.10*(TOffLocationsPromPCA - (windowStart + windoWBeginSample)));
    [~,temp1] = min(neonateECGPromPCA(windowStart + windoWBeginSample:TOffLocationsPromPCA - diffValue));
    TPeakLocationsPromPCA = windowStart + windoWBeginSample + temp1 - 1;

    TPeakValuePromPCA = neonateECGPromPCA(TPeakLocationsPromPCA);            
    TOffValuePromPCA  = neonateECGPromPCA(TOffLocationsPromPCA);

end


TSlope1 = mean(diff(neonateECGPromPCA(TPeakLocationsPromPCA:TOffLocationsPromPCA)));

TSlope2 = (neonateECGPromPCA(TOffLocationsPromPCA) - neonateECGPromPCA(TPeakLocationsPromPCA))/(TOffLocationsPromPCA - TPeakLocationsPromPCA);

RPeaksLocation = neonateECGPromPCAPeakLocUpdated2;
RPeaks = neonateECGPromPCA(RPeaksLocation);
TbyRRatio = TPeakValuePromPCA/RPeaks;

normQTOffDuration = (TOffLocationsPromPCA - RPeaksLocation)/sampFreq;
normQTDuration = (TPeakLocationsPromPCA - RPeaksLocation)/sampFreq;
normTTOffDuration = (TOffLocationsPromPCA - TPeakLocationsPromPCA)/sampFreq;
normEnergyPromPCA = sum(abs(neonateECGPromPCA).^2);
normEnergyQTOffPromPCA = sum(abs(neonateECGPromPCA(RPeaksLocation:TOffLocationsPromPCA)).^2);

TFeatures.TbyRRatio = TbyRRatio;
TFeatures.TPeakValuePromPCA = TPeakValuePromPCA;
TFeatures.TPeakLocationsPromPCA = TPeakLocationsPromPCA;
TFeatures.TOffValuePromPCA = TOffValuePromPCA;
TFeatures.TWaveShapePromPCA = TWaveShapePromPCA;
TFeatures.TOffLocationsPromPCA = TOffLocationsPromPCA;
TFeatures.TSlope1 = TSlope1;
TFeatures.TSlope2 = TSlope2;


TFeatures.normQTOffDuration = normQTOffDuration;
TFeatures.normQTDuration = normQTDuration;
TFeatures.normTTOffDuration = normTTOffDuration;
TFeatures.normEnergyPromPCA = normEnergyPromPCA;
TFeatures.normEnergyQTOffPromPCA = normEnergyQTOffPromPCA;




