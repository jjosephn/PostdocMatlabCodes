function [signalRFeatures] = panTompkinNeo(inSignal,sampFreq)

% % %  Pan-Tompkins Algorithm. Customized for neonatal ECG.

inSignal = inSignal/max(abs(inSignal));

% %% Low Pass Filter
% bLowPass = zeros(1,13);
% bLowPass(1) = 1;
% bLowPass(7) = -2;
% bLowPass(13) = 1;
% aLowPass = [1 -2 1];
% [signalLowPass,zfL] = filter(bLowPass,aLowPass,inSignal);
% 
% % signalLowPass = signalLowPass - mean(signalLowPass);
% signalLowPass = signalLowPass/max(abs(signalLowPass));
% 
% %% High Pass Filter
% bHighPass = zeros(1,33);
% bHighPass(1) = -1;
% bHighPass(17) = 32;
% bHighPass(33) = 1;
% aHighPass = [1 1];
% [signalHighPass,zfH] = filter(bHighPass,aHighPass,signalLowPass);
% 
% % signalHighPass = signalHighPass - mean(signalHighPass);
% signalHighPass = signalHighPass/max(abs(signalHighPass));
% 
% signalBandPass = signalHighPass;

%% Band Pass filter

cutOffFreqLow = 5;
cutOffFreqHigh = 15;
butterOrder = 1;
normCutOffFreq = [cutOffFreqLow cutOffFreqHigh]*2/sampFreq;                % between [0,1]
[bCoeff,aCoeff]= butter(butterOrder,normCutOffFreq);
signalBandPass = filtfilt(bCoeff,aCoeff,inSignal);                         % estimated baseline
signalBandPass = signalBandPass/max(abs(signalBandPass));

%% Derivative
bDerivative = [1 2 0 -2 1]*sampFreq/8;
aDerivative = 1;
signalDerivative = filtfilt(bDerivative,aDerivative,signalBandPass);

signalDerivative = signalDerivative/max(abs(signalDerivative));

%% Squaring
signalSquare = signalDerivative.^2;

%% Moving Window Integration (MVI)
QRSInitialWidth = 0.080;                                                   % Initial QRS width in ms for fetal ECG (80 -120 ms ???).
windowLengthQRS = round(QRSInitialWidth*sampFreq);
if mod(windowLengthQRS,2)==0
    windowLengthQRS = windowLengthQRS + 1;
end

signalMVI = conv(signalSquare,ones(1,windowLengthQRS)/windowLengthQRS);

signalMVI = signalMVI/max(abs(signalMVI));

%% Detect peaks
RRInitialWidth = 0.300;                                                    % Initial R-R width in ms for fetal ECG (min 300ms ~= 200 bpm ???).
minRRInterval = round(RRInitialWidth*sampFreq);
[peaksSignalMVITemp,peaksLocationMVITemp] = findpeaks(signalMVI,'minpeakdistance',minRRInterval);


%% Learning Phase 1 (LP1)
LPWidth1 = 2;                                                              % Two seconds for learning phase 1
signalThresholdMVI1 = max(signalMVI(1:LPWidth1*sampFreq))*1/3;             % Initializing 1st thresholds from two seconds of signal (THRESHOLD I1)
signalThresholdBandPass1 = max(signalBandPass(1:LPWidth1*sampFreq))*1/3;   % THRESHOLD F1

signalThresholdMVI2 = signalThresholdMVI1*1/2;                             % Initializing 2nd threshold from threshold 1 (0.5*threshold 2)   (THRESHOLD I2)  
signalThresholdBandPass2 = signalThresholdBandPass1*1/2;                   % THRESHOLD F2

signalPeakMVI = signalThresholdMVI1;                                       % Initial values of signal peak MVI (SPKI)
signalPeakBandPass = signalThresholdBandPass1;                             % Initial values of signal peak Band Pasd (SPKF) 

noiseThresholdMVI = mean(signalMVI(1:LPWidth1*sampFreq))*2/3;
noiseThresholdBandPass = mean(signalBandPass(1:LPWidth1*sampFreq))*2/3;

noisePeakMVI = noiseThresholdMVI;                                          % Initial values of noise peaks
noisePeakBandPass = noiseThresholdBandPass;

peaksSignalBandPassTemp = zeros(1,length(peaksSignalMVITemp));
peaksLocationBandPassTemp = zeros(1,length(peaksSignalMVITemp));

peaksSignalBandPass = zeros(1,length(peaksSignalMVITemp));
peaksLocationBandPass = zeros(1,length(peaksSignalMVITemp));

peaksSignalMVI = zeros(1,length(peaksSignalMVITemp));
peaksLocationMVI = zeros(1,length(peaksSignalMVITemp));

detectedPeakMVI = zeros(1,length(peaksSignalMVITemp));
detectedPeakBandPass = zeros(1,length(peaksSignalMVITemp));

QRSComplexNum = 0;
bufferSize = 8;
selectedRRBandPass = zeros(1,bufferSize);
selectedRRBandPassIndex = 1;
selectedRRBandPassFlag = 0;
QRSSlopeWidth = QRSInitialWidth;
QRSSlopeWidthSamples = round(QRSSlopeWidth*sampFreq);
TWaveEnd = 0.360;
peakNos = 1;
searchFlagOn = 1;
currentQRSComplexNum = 0;                            % Locating peaks in bandpass filtered signal
while peakNos <= length(peaksSignalMVITemp)  
    if peaksLocationMVITemp(peakNos) - windowLengthQRS > 0 && peaksLocationMVITemp(peakNos) <= length(signalBandPass) - windowLengthQRS
        [peaksSignalBandPassTemp(peakNos),peaksLocationBandPassTemp(peakNos)] = max(signalBandPass(peaksLocationMVITemp(peakNos) - windowLengthQRS ... 
        :peaksLocationMVITemp(peakNos) + windowLengthQRS));
        peaksLocationBandPassTemp(peakNos) = peaksLocationMVITemp(peakNos) - windowLengthQRS - 1 + peaksLocationBandPassTemp(peakNos);
    elseif peaksLocationMVITemp(peakNos) - windowLengthQRS <= 0
        [peaksSignalBandPassTemp(peakNos),peaksLocationBandPassTemp(peakNos)] = max(signalBandPass(1:peaksLocationMVITemp(peakNos) + windowLengthQRS));
        peaksLocationBandPassTemp(peakNos) = peaksLocationBandPassTemp(peakNos);
    elseif peaksLocationMVITemp(peakNos) + windowLengthQRS > length(signalBandPass)
        [peaksSignalBandPassTemp(peakNos),peaksLocationBandPassTemp(peakNos)] = max(signalBandPass(peaksLocationMVITemp(peakNos) - windowLengthQRS:end));
        peaksLocationBandPassTemp(peakNos) = peaksLocationMVITemp(peakNos) - windowLengthQRS - 1  + peaksLocationBandPassTemp(peakNos);
    end
    
    % Apply thresholds for MVI signal
    if peaksSignalMVITemp(peakNos) > signalThresholdMVI1
        detectedPeakMVI(peakNos) = 1;
        signalPeakMVI = 0.125*peaksSignalMVITemp(peakNos) + 0.875*signalPeakMVI;
    elseif peaksSignalMVITemp(peakNos) > signalThresholdMVI2 && peaksSignalMVITemp(peakNos) <= signalThresholdMVI1
        detectedPeakMVI(peakNos) = 1;                                      % Check for peaks using second threshold
        signalPeakMVI = 0.25*peaksSignalMVITemp(peakNos) + 0.75*signalPeakMVI;
    else
        detectedPeakMVI(peakNos) = 0;
        noisePeakMVI = 0.125*peaksSignalMVITemp(peakNos) + 0.875*noisePeakMVI;
    end
    
    % Apply thresholds for Band Pass signal
    if peaksSignalBandPassTemp(peakNos) > signalThresholdBandPass1
        detectedPeakBandPass(peakNos) = 1;
        signalPeakBandPass = 0.125*peaksSignalBandPassTemp(peakNos) + 0.875*signalPeakBandPass;
    elseif peaksSignalBandPassTemp(peakNos) > signalThresholdBandPass2 && peaksSignalBandPassTemp(peakNos) <= signalThresholdBandPass1
        detectedPeakBandPass(peakNos) = 1;                                 % Check for peaks using second threshold
        signalPeakBandPass = 0.25*peaksSignalBandPassTemp(peakNos) + 0.75*signalPeakBandPass;        
    else
        noisePeakBandPass = 0.125*peaksSignalBandPassTemp(peakNos) + 0.875*noisePeakBandPass;
    end
    
    
    % Compare detected peaks of MVI and Band Pass signal
    if detectedPeakMVI(peakNos) && detectedPeakBandPass(peakNos)        
        QRSComplexNum = QRSComplexNum + 1;
        peaksSignalBandPass(QRSComplexNum) = peaksSignalBandPassTemp(peakNos);
        peaksLocationBandPass(QRSComplexNum) = peaksLocationBandPassTemp(peakNos);

        peaksSignalMVI(QRSComplexNum) = peaksSignalMVITemp(peakNos);
        peaksLocationMVI(QRSComplexNum) = peaksLocationMVITemp(peakNos);
        
        if QRSComplexNum == 3
            RRAverageBandPass2 = mean(diff(peaksLocationBandPass(1:3))/sampFreq);
            
            RRLowLimitBandPass = 0.92*RRAverageBandPass2;
            RRHighLimitBandPass = 1.16*RRAverageBandPass2;
            RRMissedLimitBandPass = 1.66*RRAverageBandPass2;
            
        elseif QRSComplexNum >= bufferSize + 1
            recentRRBandPass = diff(peaksLocationBandPass(QRSComplexNum - bufferSize:QRSComplexNum))/sampFreq;        
            RRAverageBandPass1 = mean(recentRRBandPass);
            currentRRBandPass = diff(peaksLocationBandPass(QRSComplexNum - 1:QRSComplexNum))/sampFreq;
            
            % Selecting RR values within low and high limits until minimum
            % requirement is achieved (= bufferSize)
            if selectedRRBandPassIndex <= bufferSize + 1                   
                if currentRRBandPass >= RRLowLimitBandPass && currentRRBandPass <= RRHighLimitBandPass
                    selectedRRBandPass(selectedRRBandPassIndex) = currentRRBandPass;
                    if selectedRRBandPassIndex == bufferSize + 1
                        selectedRRBandPassFlag = 1;
                    end
                    selectedRRBandPassIndex = selectedRRBandPassIndex + 1;
                end
            end
            
            % Selecting bufferSize number of RR values within the limits and compute RR Average 2 
            if currentRRBandPass >= RRLowLimitBandPass && currentRRBandPass <= RRHighLimitBandPass && selectedRRBandPassFlag
                selectedRRBandPass = [selectedRRBandPass(2:bufferSize) currentRRBandPass];
                RRAverageBandPass2 = mean(selectedRRBandPass);
            end
            
            RRLowLimitBandPass = 0.92*RRAverageBandPass2;
            RRHighLimitBandPass = 1.16*RRAverageBandPass2;
            RRMissedLimitBandPass = 1.66*RRAverageBandPass2;
            
        end
    end    
    
    if QRSComplexNum >= 3
        TWaveEndLocation = round(TWaveEnd*sampFreq);
        currentRRBandPass = (peaksLocationBandPass(QRSComplexNum) - peaksLocationBandPass(QRSComplexNum - 1))/sampFreq;
        previousRRBandPass = (peaksLocationBandPass(QRSComplexNum - 1) - peaksLocationBandPass(QRSComplexNum - 2))/sampFreq;
        
        currentPeakLocationBandPass = peaksLocationBandPass(QRSComplexNum);
        previousPeakLocationBandPass = peaksLocationBandPass(QRSComplexNum - 1);
        if currentPeakLocationBandPass - previousPeakLocationBandPass <= TWaveEndLocation
            currentSlopeBandPass = mean(diff(signalBandPass(peaksLocationBandPass(QRSComplexNum) - QRSSlopeWidthSamples:peaksLocationBandPass(QRSComplexNum))));      % mean slope of current R wave
            previousSlopeBandPass = mean(diff(signalBandPass(peaksLocationBandPass(QRSComplexNum - 1) - QRSSlopeWidthSamples:peaksLocationBandPass(QRSComplexNum - 1)))); % mean slope of previous R wave
            if currentSlopeBandPass <= 0.5*previousSlopeBandPass           % Explore possible T wave.
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                currentRRBandPass = (peaksLocationBandPass(QRSComplexNum) - peaksLocationBandPass(QRSComplexNum - 2))/sampFreq;
                previousRRBandPass = (peaksLocationBandPass(QRSComplexNum - 1) - peaksLocationBandPass(QRSComplexNum - 2))/sampFreq;
                
                peaksSignalBandPass(QRSComplexNum) = [];
                peaksLocationBandPass(QRSComplexNum) = [];

                peaksSignalMVI(QRSComplexNum) = [];
                peaksLocationMVI(QRSComplexNum) = [];
                
                QRSComplexNum = QRSComplexNum - 1;
                
                noisePeakMVI = 0.125*peaksSignalMVITemp(peakNos) + 0.875*noisePeakMVI;
                noisePeakBandPass = 0.125*peaksSignalBandPassTemp(peakNos) + 0.875*noisePeakBandPass;                
            end
        end
        
        if abs(currentRRBandPass - previousRRBandPass) >= RRMissedLimitBandPass           % Search back for missed R peak
           [searchBackPeakBandPassTemp,searchBackPeaksLocationBandPassTemp] = max(signalBandPass(peaksLocationBandPass(QRSComplexNum - 1):peaksLocationBandPass(QRSComplexNum)));
           searchBackPeaksLocationBandPassTemp = searchBackPeaksLocationBandPassTemp + peaksLocationBandPass(QRSComplexNum - 1);
           if searchBackPeakBandPassTemp > signalThresholdBandPass2 && searchBackPeakBandPassTemp <= signalThresholdBandPass1
               peaksSignalBandPass(QRSComplexNum) = searchBackPeakBandPassTemp;
               peaksLocationBandPass(QRSComplexNum) = searchBackPeaksLocationBandPassTemp;
               %%%%%%%%%%%%
               signalPeakBandPass = 0.25*peaksSignalBandPassTemp(peakNos) + 0.75*signalPeakBandPass;
           end               
        end
    end
  

    if QRSComplexNum > 3   
        twoRRMean = mean(diff(peaksLocationBandPass(QRSComplexNum - 3 + 1:QRSComplexNum)));
        if peaksLocationMVITemp(peakNos) - peaksLocationBandPass(QRSComplexNum) > twoRRMean && searchFlagOn && (peaksLocationMVITemp(peakNos) + LPWidth1*sampFreq) < length(inSignal)
            currentQRSComplexNum = QRSComplexNum;
            tempPeakLocationMVI = find(peaksLocationMVITemp>peaksLocationBandPass(QRSComplexNum) + round(twoRRMean)/2);
            peakNos = tempPeakLocationMVI(1);
            signalThresholdMVI1 = max(signalMVI(peaksLocationMVITemp(peakNos):peaksLocationMVITemp(peakNos) + LPWidth1*sampFreq))*1/3;             % Initializing 1st thresholds from two seconds of signal (THRESHOLD I1)
            signalThresholdBandPass1 = max(signalBandPass(peaksLocationMVITemp(peakNos):peaksLocationMVITemp(peakNos) + LPWidth1*sampFreq))*1/3;   % THRESHOLD F1

            signalThresholdMVI2 = signalThresholdMVI1*1/2;                             % Initializing 2nd threshold from threshold 1 (0.5*threshold 2)   (THRESHOLD I2)  
            signalThresholdBandPass2 = signalThresholdBandPass1*1/2;                   % THRESHOLD F2        
            
            signalPeakMVI = signalThresholdMVI1;                                       % Initial values of signal peak MVI (SPKI)
            signalPeakBandPass = signalThresholdBandPass1;                             % Initial values of signal peak Band Pasd (SPKF) 

            noiseThresholdMVI = mean(signalMVI(peaksLocationMVITemp(peakNos):peaksLocationMVITemp(peakNos) + LPWidth1*sampFreq))*2/3;
            noiseThresholdBandPass = mean(signalBandPass(peaksLocationMVITemp(peakNos):peaksLocationMVITemp(peakNos) + LPWidth1*sampFreq))*2/3;

            noisePeakMVI = noiseThresholdMVI;                                          % Initial values of noise peaks
            noisePeakBandPass = noiseThresholdBandPass;
            searchFlagOn = 0;

        else
            signalThresholdMVI1 = noisePeakMVI + 0.25*(signalPeakMVI - noisePeakMVI);
            signalThresholdMVI2 = signalThresholdMVI1*1/2;

            signalThresholdBandPass1 = noisePeakBandPass + 0.25*(signalPeakBandPass - noisePeakBandPass);
            signalThresholdBandPass2 = signalThresholdBandPass1*1/2;    
            peakNos = peakNos + 1;
            if QRSComplexNum > currentQRSComplexNum
                searchFlagOn = 1;
            end
        end
    else
        signalThresholdMVI1 = noisePeakMVI + 0.25*(signalPeakMVI - noisePeakMVI);
        signalThresholdMVI2 = signalThresholdMVI1*1/2;

        signalThresholdBandPass1 = noisePeakBandPass + 0.25*(signalPeakBandPass - noisePeakBandPass);
        signalThresholdBandPass2 = signalThresholdBandPass1*1/2;
        peakNos = peakNos + 1;
        if QRSComplexNum > currentQRSComplexNum
            searchFlagOn = 1;
        end
    end
  
    
end

peaksSignalBandPass(peaksSignalBandPass==0) = [];
peaksLocationBandPass(peaksLocationBandPass==0) = [];

peaksSignalMVI(peaksSignalMVI==0) = [];
peaksLocationMVI(peaksLocationMVI==0) = [];

peaksLocationBandPass = peaksLocationBandPass(peaksSignalBandPass > 0.5*mean(peaksSignalBandPass));
peaksSignalBandPass = peaksSignalBandPass(peaksSignalBandPass > 0.5*mean(peaksSignalBandPass));


%% Correction

peaksLocationSignal = peaksLocationBandPass;
peaksValueSignal = zeros(1,length(peaksLocationSignal));

if peaksLocationBandPass(1) - 3 <= 0
    [peaksValueSignal(1),tempLoc] = max(inSignal(1:peaksLocationSignal(1) + 3));
    peaksLocationSignal(1) = tempLoc;
else
    [peaksValueSignal(1),tempLoc] = max(inSignal(peaksLocationSignal(1) - 3:peaksLocationSignal(1) + 3));
    peaksLocationSignal(1) = peaksLocationSignal(1) - 3 - 1 + tempLoc;
end


for var1=2:1:length(peaksLocationSignal) - 1
    [peaksValueSignal(var1),tempLoc] = max(inSignal(peaksLocationSignal(var1) - 3:peaksLocationSignal(var1) + 3));
    peaksLocationSignal(var1) = peaksLocationSignal(var1) - 3 - 1 + tempLoc;
end

if peaksLocationBandPass(end) + 3 > length(inSignal)
    [peaksValueSignal(end),tempLoc] = max(inSignal(peaksLocationSignal(end) - 3:end));
    peaksLocationSignal(end) = peaksLocationSignal(end) - 3 - 1 + tempLoc;
else    
    [peaksValueSignal(end),tempLoc] = max(inSignal(peaksLocationSignal(end) - 3:peaksLocationSignal(end) + 3));
    peaksLocationSignal(end) = peaksLocationSignal(end) - 3 - 1 + tempLoc;
end


signalRFeatures.signalBandPass = signalBandPass;
signalRFeatures.peaksSignalBandPass = peaksSignalBandPass;
signalRFeatures.peaksLocationBandPass = peaksLocationBandPass;

signalRFeatures.signal = inSignal;
signalRFeatures.peaksValueSignal = peaksValueSignal;
signalRFeatures.peaksLocationSignal = peaksLocationSignal;

signalRFeatures.peaksSignalMVI = peaksSignalMVI;
signalRFeatures.peaksLocationMVI = peaksLocationMVI;

signalRFeatures.RRAverageBandPass1 = RRAverageBandPass1;
signalRFeatures.RRAverageBandPass2 = RRAverageBandPass2;



