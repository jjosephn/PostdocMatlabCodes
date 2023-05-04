function [TFeatures] = TFeaturesNeo(inSignal,QRSFeatures,sampFreq)



TFeatures.inSignal = inSignal;

signalBandPass = QRSFeatures.signalBandPass;
peaksSignalBandPass = QRSFeatures.peaksSignalBandPass;
peaksLocationBandPass = QRSFeatures.peaksLocationBandPass;

inSignal = signalBandPass;
inSignal = inSignal - mean(inSignal);
inSignal = inSignal/max(abs(inSignal));

QRSInitialWidth = 0.080;                                                   % Initial QRS width in ms for fetal ECG (80 -120 ms ???).
windowLengthQRS = round(QRSInitialWidth*sampFreq);
if mod(windowLengthQRS,2) == 1
    windowLengthQRS = windowLengthQRS + 1;
end

paddedWindowLength = windowLengthQRS;
paddedInSignal = [zeros(1,paddedWindowLength) inSignal];

%% T wave estimation


% Low pass filter

cutOffFreqLowTWave = 50; %
butterOrderTWave = 1;
normCutOffFreqTWave = cutOffFreqLowTWave/(sampFreq/2); % between [0,1]
[bCoeffTWave,aCoeffTWave]= butter(butterOrderTWave,normCutOffFreqTWave,'low');
signalBandPassTWave = filtfilt(bCoeffTWave,aCoeffTWave,paddedInSignal);                         % estimated baseline
signalBandPassTWave = signalBandPassTWave(paddedWindowLength + 1:end);

signalBandPassTWave = signalBandPassTWave - mean(signalBandPassTWave);
signalBandPassTWave = signalBandPassTWave/max(abs(signalBandPassTWave));

smoothSignalBandPassTWave = smooth(signalBandPassTWave,windowLengthQRS + 1,'sgolay',5);
smoothSignalBandPassTWave = smoothSignalBandPassTWave';

smoothSignalBandPassTWave = smoothSignalBandPassTWave - mean(smoothSignalBandPassTWave);
smoothSignalBandPassTWave = smoothSignalBandPassTWave/max(abs(smoothSignalBandPassTWave));

derivativeTWaveSignal = diff(smoothSignalBandPassTWave);
derivativeTWaveSignal = [derivativeTWaveSignal 0];

meanRRInterval = mean(diff(peaksLocationBandPass)/sampFreq);

derivativeTWaveSignal = derivativeTWaveSignal - mean(derivativeTWaveSignal);
derivativeTWaveSignal = derivativeTWaveSignal/max(abs(derivativeTWaveSignal));
windoWBeginSample = round(0.2*meanRRInterval*sampFreq);
windowMiddleSample1 = round(0.25*meanRRInterval*sampFreq);
windowMiddleSample2 = round(0.50*meanRRInterval*sampFreq);
windowEndSample = round(0.65*meanRRInterval*sampFreq);


TWaveShape = zeros(1,length(peaksLocationBandPass));
TPeakValue = zeros(1,length(peaksLocationBandPass));
TPeakLocations = zeros(1,length(peaksLocationBandPass));
TOffValue = zeros(1,length(peaksLocationBandPass));
TOffLocations = zeros(1,length(peaksLocationBandPass));
QTPeaks = zeros(2,length(peaksLocationBandPass));
QTLocations = zeros(2,length(peaksLocationBandPass));

currentWindow = zeros(1,length(derivativeTWaveSignal));
for var1=1:1:length(peaksLocationBandPass)
    windowStart = peaksLocationBandPass(var1);    
    if windowStart + windowEndSample <= length(derivativeTWaveSignal)
                
        for var2 = windowStart + windoWBeginSample:1:windowStart + windowMiddleSample1
            currentWindow(var2) = ((var2  - windowStart - windoWBeginSample)/(windowMiddleSample1 - windoWBeginSample));
        end
        currentWindow(windowStart + windowMiddleSample1 + 1:windowStart + windowMiddleSample2) = 1;
        for var2 = windowStart + windowMiddleSample2 + 1:1:windowStart + windowEndSample
            currentWindow(var2) = 1 - ((var2 - windowStart - windowMiddleSample2)/(windowEndSample - windowMiddleSample2));
        end
        
        windowWidth = find(currentWindow>0);
        winDerivativeTWaveSignal = derivativeTWaveSignal.*currentWindow;
        
        [maxValue,maxValueLocation] = max(winDerivativeTWaveSignal(windowStart + windoWBeginSample:windowStart + windowEndSample));
        [minValue,minValueLocation] = min(winDerivativeTWaveSignal(windowStart + windoWBeginSample:windowStart + windowEndSample));
        
        maxValueLocation = windowStart + windoWBeginSample + maxValueLocation - 1;
        minValueLocation = windowStart + windoWBeginSample + minValueLocation - 1;

        if maxValueLocation < minValueLocation
            TWaveShape(var1) = 1;
            thresholdTWave = minValue/2;
            
            var3 = minValueLocation;
            var4 = 1;
            thresholdFlag = 1;
            while thresholdFlag
                currentValue = winDerivativeTWaveSignal(var3);
                if var3 == windowWidth(length(windowWidth))
                    TOffLocations(var1) = minValueLocation + var4 - 1;
                    thresholdFlag = 0;
                elseif currentValue > thresholdTWave
                    TOffLocations(var1) = minValueLocation + var4 - 1;
                    thresholdFlag = 0;
                end
                var3 = var3 + 1;
                var4 = var4 + 1;
            end
            
            diffValue = round(0.10*(TOffLocations(var1) - (windowStart + windoWBeginSample)));
            [~,temp1] = max(inSignal(windowStart + windoWBeginSample:TOffLocations(var1) - diffValue));
            TPeakLocations(var1) = windowStart + windoWBeginSample + temp1 - 1;
            
            TPeakValue(var1) = inSignal(TPeakLocations(var1));            
            TOffValue(var1)  = inSignal(TOffLocations(var1));            
            QTPeaks(1,var1) = peaksSignalBandPass(var1);
            QTPeaks(2,var1) = TPeakValue(var1);
            QTLocations(1,var1) = peaksLocationBandPass(var1);
            QTLocations(2,var1) = TPeakLocations(var1);               
            
        elseif maxValueLocation > minValueLocation        
            TWaveShape(var1) = -1;
            thresholdTWave = maxValue/2;
            
            var3 = maxValueLocation;
            var4 = 1;
            thresholdFlag = 1;
            while thresholdFlag
                currentValue = winDerivativeTWaveSignal(var3);
                if var3 == windowWidth(length(windowWidth))
                    TOffLocations(var1) = maxValueLocation + var4 - 1;
                    thresholdFlag = 0;
                elseif currentValue < thresholdTWave
                    TOffLocations(var1) = maxValueLocation + var4 - 1;
                    thresholdFlag = 0;
                end
                var3 = var3 + 1;
                var4 = var4 + 1;
            end
            
            diffValue = round(0.10*(TOffLocations(var1) - (windowStart + windoWBeginSample)));
            [~,temp1] = min(inSignal(windowStart + windoWBeginSample:TOffLocations(var1) - diffValue));
            TPeakLocations(var1) = windowStart + windoWBeginSample + temp1 - 1;
            
            TPeakValue(var1) = inSignal(TPeakLocations(var1));            
            TOffValue(var1)  = inSignal(TOffLocations(var1));
            QTPeaks(1,var1) = peaksSignalBandPass(var1);
            QTPeaks(2,var1) = TPeakValue(var1);
            QTLocations(1,var1) = peaksLocationBandPass(var1);
            QTLocations(2,var1) = TPeakLocations(var1);               
            
        end
    end    
    currentWindow = zeros(1,length(derivativeTWaveSignal));
end

QTPeaks(:,TPeakLocations==0) = [];
QTLocations(:,TPeakLocations==0) = [];    

TPeakValue(TPeakLocations==0) = [];
TPeakLocations(TPeakLocations==0) = [];
TOffValue(TOffLocations==0) = [];
TWaveShape(TOffLocations==0) = [];
TOffLocations(TOffLocations==0) = [];

TFeatures.TPeakValue = TPeakValue;
TFeatures.TPeakLocations = TPeakLocations;
TFeatures.TOffValue = TOffValue;
TFeatures.TWaveShape = TWaveShape;
TFeatures.TOffLocations = TOffLocations;
TFeatures.QTPeaks = QTPeaks;
TFeatures.QTLocations = QTLocations;



