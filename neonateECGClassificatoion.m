clc;
close all;
clearvars; 

load 'neonateECGSegment600Sec';

%%%%%%% good samples

includeSet1HighpH1 = [2 5 35 36 43 44 79 81 86 89 90 91 92 93 94 96 97 98 100 ... 
    111 131 210 211 212 213 215 216 217 218 219 226 227 230 231 232 235 ... 
    236 239 240 256];                                          

includeSet1PosTHighpH1 = [1 40 56 57 58 107 199 200 201 202 248];                    % positive T-wave

includeSet1HighpH2 = [23 25 26 27 28 30 31 41 42 45 46 47 48 50 ...
    51 109 110 119 122 127 142 156 157 159 164 ...
    165 167 171 193 194 195 251];

includeSet1PosTHighpH2 = [13 14 15 16 71 72 73 77 78 117 125 149 158 161 163 187 188];                                % positive T-wave

includeSet1LowpH2 = [162 166 174 175];

includeSet1PosTLowpH2 = 76;                                                 % positive T-wave

%%%%%% Less than 60 seconds

includeSet2HighpH1 = [3 37 53 82 83 84 88 95 99 214 262 270];

includeSet2PosTHighpH1 = 198;                                               % positive T-wave

includeSet2HighpH2 = 29;

includeSet2LowpH2 = [];

includeSet3LowpH1 = [145 146];

includeSet3LowpH2 = 173;

includeSet3PosTLowpH2 = 176;

positiveSet = [includeSet1PosTHighpH1 includeSet1PosTHighpH2 includeSet1PosTLowpH2 includeSet2PosTHighpH1 includeSet3PosTLowpH2];

% %%%%%%% Whole data

% includeSet1 = [includeSet1HighpH1 includeSet1PosTHighpH1 includeSet2HighpH1 includeSet2PosTHighpH1 includeSet3LowpH1];
% includeSet2 = [includeSet1HighpH2 includeSet1PosTHighpH2 includeSet2HighpH2 includeSet1LowpH2 includeSet1PosTLowpH2 includeSet2LowpH2 includeSet3LowpH2 includeSet3PosTLowpH2];

%%%%%%% Whole data - positive T-wave

% includeSet1 = [includeSet1HighpH1 includeSet2HighpH1 includeSet3LowpH1];
% includeSet2 = [includeSet1HighpH2 includeSet2HighpH2 includeSet1LowpH2 includeSet2LowpH2 includeSet3LowpH2];

%%%%%%% Low pH 

% includeSet1 = [includeSet3LowpH1];
% includeSet2 = [includeSet1LowpH2 includeSet1PosTLowpH2 includeSet2LowpH2 includeSet3LowpH2 includeSet3PosTLowpH2];

%%%%%%% Low pH  - positive T-wave

% includeSet1 = [includeSet3LowpH1];
% includeSet2 = [includeSet1LowpH2 includeSet2LowpH2 includeSet3LowpH2];

%%%%%%% High pH 

% includeSet1 = [includeSet1HighpH1 includeSet1PosTHighpH1 includeSet2HighpH1 includeSet2PosTHighpH1];
% includeSet2 = [includeSet1HighpH2 includeSet1PosTHighpH2 includeSet2HighpH2];

%%%%%%% High pH  - positive T-wave 

% includeSet1 = [includeSet1HighpH1 includeSet2HighpH1];
% includeSet2 = [includeSet1HighpH2 includeSet2HighpH2];

%%%%%%%  Equal Sets

% includeSet1 = [145 146 ... 
%     2 35 94 216];
% includeSet2 = [162 166 173 174 175 176 ... 
%     26 27 167 195];

% %%%%%%  Equal Sets 2
% 
% includeSet1 = [145 146 ... 
%     2 5 51 89 91 131];
% includeSet2 = [162 166 173 174 175 176 ... 
%     45 149];

%%%%%%%  Equal Sets random 2

% includeSet1 = [145 146 ... 
%     2 91 131 201 230 240 ...
%     5 89 97 100 248];
% includeSet2 = [76 162 166 173 174 175 176 ... 
%     149 163 171 ...
%     45 51 119];

%%%%%%%  Equal Sets random - pos T

includeSet1 = [145 146 ... 
    2 35 43 94 230];
includeSet2 = [162 166 173 174 175 ... 
    27 40];

%%%%%%%  Equal Sets random

% includeSet1 = [145 146 ... 
%     2 35 90 94 230];
% includeSet2 = [76 162 166 173 174 175 176 ... 
%     27 40 72 149];

excludeSet = [1 14 15 16 25 77 107 109 110 122 157 159 210 262 ... 
    73];

factor = 1;
% while factor <= 5
segmentNameAll = fieldnames(neonateECGSegment);
meanTbyQRSRatio = zeros(1,length(fieldnames(neonateECGSegment)));
meanTSlope = zeros(1,length(fieldnames(neonateECGSegment)));
segmentNoArray = zeros(1,length(fieldnames(neonateECGSegment)));
pH = zeros(1,length(fieldnames(neonateECGSegment)));
TShapePercentage = zeros(1,length(fieldnames(neonateECGSegment)));
RRAvg = zeros(1,length(fieldnames(neonateECGSegment)));
heartRate = zeros(1,length(fieldnames(neonateECGSegment)));
meanT = zeros(1,length(fieldnames(neonateECGSegment)));
medianT = zeros(1,length(fieldnames(neonateECGSegment)));
varianceT = zeros(1,length(fieldnames(neonateECGSegment)));

normQTOffDuration = zeros(1,length(fieldnames(neonateECGSegment)));
normQTDuration = zeros(1,length(fieldnames(neonateECGSegment)));
normTTOffDuration = zeros(1,length(fieldnames(neonateECGSegment)));
normEnergyPromPCA = zeros(1,length(fieldnames(neonateECGSegment)));
normEnergyQTOffPromPCA = zeros(1,length(fieldnames(neonateECGSegment)));


for segmentNum=1:1:length(fieldnames(neonateECGSegment))
    fprintf('%s \n',segmentNameAll{segmentNum});
    sampFreq = neonateECGSegment.(segmentNameAll{segmentNum}).sampFreq;
    segmentMinutes = 1/6;
    segmentLength = segmentMinutes*60*sampFreq;

    segmentSamplePoint = neonateECGSegment.(segmentNameAll{segmentNum}).segmentSamplePoint;
    
    start = segmentSamplePoint - factor*segmentLength + 1;
    stop = segmentSamplePoint - (factor - 1)*segmentLength;

    if (~isempty(includeSet1(includeSet1==segmentNum)) || ~isempty(includeSet2(includeSet2==segmentNum))) && isempty(excludeSet(excludeSet==segmentNum))
        if length(neonateECGSegment.(segmentNameAll{segmentNum}).record(:,1)) > 1
            if ~isempty(includeSet1(includeSet1==segmentNum))
                neonateECG = neonateECGSegment.(segmentNameAll{segmentNum}).record(1,start:stop);
            elseif ~isempty(includeSet2(includeSet2==segmentNum))
                neonateECG = neonateECGSegment.(segmentNameAll{segmentNum}).record(2,start:stop);
            end
        else
            neonateECG = neonateECGSegment.(segmentNameAll{segmentNum}).record(1,start:stop);
        end


        QRSInitialWidth = 0.080;                                                   % Initial QRS width in ms for fetal ECG (80 -120 ms ???).
        windowLengthQRS = round(QRSInitialWidth*sampFreq);
        if mod(windowLengthQRS,2) == 1
            windowLengthQRS = windowLengthQRS + 1;
        end

        paddedWindowLength = windowLengthQRS;
        paddedInSignal = [zeros(1,paddedWindowLength) neonateECG zeros(1,paddedWindowLength)];

        neonateECG = paddedInSignal;

        signalBLRemoved = baselineRemoval(neonateECG,sampFreq);
        signalNoiseRemoved = highFreqNoiseRemoval(signalBLRemoved,sampFreq,1);

        signalNoiseRemoved = signalNoiseRemoved(paddedWindowLength + 1:length(paddedInSignal) - paddedWindowLength);

        preProcessedSignal = smooth(signalNoiseRemoved,windowLengthQRS + 1,'sgolay',5);
        preProcessedSignal = preProcessedSignal';

        preProcessedSignal(1,:) = preProcessedSignal(1,:)/max(abs(preProcessedSignal(1,:)));
        preProcessedSignal(1,:) = preProcessedSignal(1,:) - mean(preProcessedSignal(1,:));     

        QRSFeatures = QRSFeaturesNeo(preProcessedSignal,sampFreq);
        TFeaturesPCA = TFeaturesPCANeo(QRSFeatures,sampFreq);

        hrSampleLen = 2;
        contHeartRate = zeros(1,length(QRSFeatures.peaksLocationBandPass));
        for var1=hrSampleLen:1:length(QRSFeatures.peaksLocationBandPass)
            contHeartRate(var1) = 60/mean(diff(QRSFeatures.peaksLocationBandPass(var1 - hrSampleLen + 1:var1))/sampFreq);
        end

        contHeartRate = contHeartRate(hrSampleLen:end);

        TPeak = TFeaturesPCA.TPeakValuePromPCA;

        TbyQRSRatio = TFeaturesPCA.TbyRRatio;
        TSlope = TFeaturesPCA.TSlope1;
        
        
        normQTOffDuration(segmentNum) = TFeaturesPCA.normQTOffDuration;
        normQTDuration(segmentNum) = TFeaturesPCA.normQTDuration;
        normTTOffDuration(segmentNum) = TFeaturesPCA.normTTOffDuration;
        normEnergyPromPCA(segmentNum) = TFeaturesPCA.normEnergyPromPCA;
        normEnergyQTOffPromPCA(segmentNum) = TFeaturesPCA.normEnergyQTOffPromPCA;
        
        segmentNoArray(segmentNum) = segmentNum;

        meanT(segmentNum) = TPeak;
        meanTbyQRSRatio(segmentNum) = TbyQRSRatio;
        meanTSlope(segmentNum) = TSlope;
        heartRate(segmentNum) = ceil(length(QRSFeatures.peaksLocationBandPass)/segmentMinutes);
        RRAvg(segmentNum) = mean(diff(QRSFeatures.peaksLocationBandPass)/sampFreq);

        pH(segmentNum) = neonateECGSegment.(segmentNameAll{segmentNum}).param.pH;
    end    
end

segmentLocation = find(pH>0);
segmentNoArray = segmentNoArray(segmentLocation);
pH = pH(segmentLocation);
meanTbyQRSRatio = meanTbyQRSRatio(segmentLocation);
meanTSlope = meanTSlope(segmentLocation);
RRAvg = RRAvg(segmentLocation);
heartRate = heartRate(segmentLocation);
meanT = meanT(segmentLocation);

normQTOffDuration = normQTOffDuration(segmentLocation);
normQTDuration = normQTDuration(segmentLocation);
normTTOffDuration = normTTOffDuration(segmentLocation);
normEnergyPromPCA = normEnergyPromPCA(segmentLocation);
normEnergyQTOffPromPCA = normEnergyQTOffPromPCA(segmentLocation);


classes = -1*ones(1,length(pH));
classes(pH<=7.2) = 1;


segmentLocation = segmentLocation';
predParam = pH;

featureVectorOriginal = [heartRate' meanTbyQRSRatio' meanTSlope' meanT'];

% featureVectorOriginal = [heartRate' meanTbyQRSRatio' meanTSlope' meanT' ... 
%     normQTOffDuration' normQTDuration' normTTOffDuration' ... 
%     normEnergyPromPCA' normEnergyQTOffPromPCA'];



weights = dec2bin(1:1:2^(length(featureVectorOriginal(1,:)))-1)-'0';

accuracy = zeros(1,length(weights(:,1)));
sensitivity = zeros(1,length(weights(:,1)));
specificity = zeros(1,length(weights(:,1)));
geometricMean = zeros(1,length(weights(:,1)));
FScore = zeros(1,length(weights(:,1)));


for featureNum = 1:1:length(weights(:,1))
    featureVector = featureVectorOriginal.*weights(featureNum,:);
    
    featureVector(:,weights(featureNum,:)==0) = [];   


    normArg = 1;
    for var2=1:1:length(featureVector(1,:))
        [featureVector(:,var2),~,~] = normalizeData(featureVector(:,var2),normArg);
    end

    inData = [predParam' featureVector];

    classesShift = classes';
    predParamShift = predParam';
    featureVectorShift = featureVector;

    rmsePredParamLin = zeros(1,length(featureVector(:,1)));
    rmsePredParamSVR = zeros(1,length(featureVector(:,1)));
    predictedSegment = zeros(1,length(featureVector(:,1)));

    truePosCount = 0;
    trueNegCount = 0;
    falsePosCount = 0;
    falseNegCount = 0;
    
    accuracyInside = zeros(1,length(featureVector(:,1)));
    sensitivityInside = zeros(1,length(featureVector(:,1)));
    specificityInside = zeros(1,length(featureVector(:,1)));
    geometricMeanInside = zeros(1,length(featureVector(:,1)));
    FScoreInside = zeros(1,length(featureVector(:,1)));

    for segmentNum = 1:1:length(featureVector(:,1))
        truePosCountInside = 0;
        trueNegCountInside = 0;
        falsePosCountInside = 0;
        falseNegCountInside = 0;


        if segmentNum==23
            fprintf('stop \n');
        end

        classesTrain = classesShift(1:end - 1);
        featureVectorTrain = featureVectorShift(1:end - 1,:);

        classesTest = classesShift(end);
        featureVectorTest = featureVectorShift(end,:);
        
%         svmmodel = svmtrain(classesTrain,featureVectorTrain,'-s 1 -t 1 -n 0.1 -q'); % Whole -T Also 0-1, 0-2, 0-3
%         svmmodel = svmtrain(classesTrain,featureVectorTrain,'-s 1 -t 1 -n 0.1 -q'); % Whole Also 0-1, 0-2, 0-3
        svmmodel = svmtrain(classesTrain,featureVectorTrain,'-s 1 -t 1 -n 0.1 -q'); % Equal - T n= 0.5 and 0.1
%         svmmodel = svmtrain(classesTrain,featureVectorTrain,'-s 1 -t 1 -n 0.1 -q'); % Equal n= 0.5 and 0.1
        [predClassSVM, ~, ~] = svmpredict(classesTest,featureVectorTest,svmmodel);

        for var2=1:1:length(predClassSVM)
            if classesTest(var2) == 1 && predClassSVM(var2) == 1
                truePosCount = truePosCount + 1;
            elseif classesTest(var2) == -1 && predClassSVM(var2) == -1
                trueNegCount = trueNegCount + 1;
            elseif classesTest(var2) == -1 && predClassSVM(var2) == 1
                falsePosCount = falsePosCount + 1;
            elseif classesTest(var2) == 1 && predClassSVM(var2) == -1
                falseNegCount = falseNegCount + 1;
            end
        end
        
        for var2=1:1:length(predClassSVM)
            if classesTest(var2) == 1 && predClassSVM(var2) == 1
                truePosCountInside = 1;
            elseif classesTest(var2) == -1 && predClassSVM(var2) == -1
                trueNegCountInside = 1;
            elseif classesTest(var2) == -1 && predClassSVM(var2) == 1
                falsePosCountInside = 1;
            elseif classesTest(var2) == 1 && predClassSVM(var2) == -1
                falseNegCountInside = 1;
            end
        end
        
        accuracyInside(segmentNum) = 100*(truePosCountInside + trueNegCountInside)/(truePosCountInside + trueNegCountInside + falsePosCountInside + falseNegCountInside);
        sensitivityInside(segmentNum) = 100*truePosCountInside/(truePosCountInside + falseNegCountInside);
        specificityInside(segmentNum) = 100*trueNegCountInside/(trueNegCountInside + falsePosCountInside);
        geometricMeanInside(segmentNum) = sqrt(sensitivityInside(segmentNum)*specificityInside(segmentNum));
        FScoreInside(segmentNum) = 100*(2*truePosCountInside)/(2*truePosCountInside + falsePosCountInside + falseNegCountInside);

        featureVectorShift = circshift(featureVectorShift,1,1);
        classesShift = circshift(classesShift,1,1);
        predParamShift = circshift(predParamShift,1,1);
        segmentLocation = circshift(segmentLocation,1,1);
    end
    
%     accuracyA(featureNum) = mean(accuracyInside);
%     sensitivity(featureNum) = mean(sensitivityInside);
%     specificity(featureNum) = mean(specificityInside);
%     geometricMean(featureNum) = mean(geometricMeanInside);
%     FScore(featureNum) = mean(FScoreInside);
    
    accuracy(featureNum) = 100*(truePosCount + trueNegCount)/(truePosCount + trueNegCount + falsePosCount + falseNegCount);
    sensitivity(featureNum) = 100*truePosCount/(truePosCount + falseNegCount);
    specificity(featureNum) = 100*trueNegCount/(trueNegCount + falsePosCount);
    geometricMean(featureNum) = sqrt(sensitivity(featureNum)*specificity(featureNum));
    FScore(featureNum) = 100*(2*truePosCount)/(2*truePosCount + falsePosCount + falseNegCount);
    fprintf(' accuracy = %2.2f %% \n sensitivity = %2.2f %% \n specificity = %2.2f %% \n geometricMean = %2.2f %% \n FScore = %2.2f %% \n', ... 
        accuracy(featureNum),sensitivity(featureNum),specificity(featureNum),geometricMean(featureNum),FScore(featureNum));
end

% barGraphData = [accuracy' sensitivity' specificity' geometricMean' FScore'];
barGraphData = [accuracy' sensitivity' specificity'];

figure
Groups = categorical({'T','TSl','TSl + T','T/QRS','T/QRS + T','T/QRS + TSl','T/QRS + TSl + T','HR', 'HR + T', ... 
    'HR + TSl','HR + TSl + T','HR + T/QRS','HR + T/QRS + T','HR + T/QRS + TSl','HR + T/QRS + TSl + T'});
Groups = reordercats(Groups,{'T','TSl','TSl + T','T/QRS','T/QRS + T','T/QRS + TSl','T/QRS + TSl + T','HR', 'HR + T', ... 
    'HR + TSl','HR + TSl + T','HR + T/QRS','HR + T/QRS + T','HR + T/QRS + TSl','HR + T/QRS + TSl + T'});
bb = bar(Groups,barGraphData);
bb(1).FaceColor = 'r';
bb(2).FaceColor = 'g';
bb(3).FaceColor = 'b';
% bb(4).FaceColor = 'm';
% bb(5).FaceColor = 'y';
xtickangle(90)
xlabel('Features')
ylabel('%')
legend({'Accuracy', 'Sensitivity', 'Specificity'},'Location','southeastoutside')
% legend({'Accuracy', 'Sensitivity', 'Specificity', 'Geometric Mean', 'F Score'},'Location','southeastoutside')
% legend('boxoff')
% ylim([0 150]);
% xlim([1 15]);

[aacVal,accInd] = max(accuracy);

fprintf('\n\n %s \naccuracy = %2.2f %% \n sensitivity = %2.2f %% \n specificity = %2.2f %% \n geometricMean = %2.2f %% \n FScore = %2.2f %% \n', ... 
        Groups(accInd), aacVal,sensitivity(accInd),specificity(accInd),geometricMean(accInd),FScore(accInd));


