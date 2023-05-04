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

% includeSet1 = includeSet3LowpH1;
% includeSet2 = [includeSet1LowpH2 includeSet1PosTLowpH2 includeSet2LowpH2 includeSet3LowpH2 includeSet3PosTLowpH2];

%%%%%%% Low pH  - positive T-wave

% includeSet1 = includeSet3LowpH1;
% includeSet2 = [includeSet1LowpH2 includeSet2LowpH2 includeSet3LowpH2];

%%%%%%% High pH 

% includeSet1 = [includeSet1HighpH1 includeSet1PosTHighpH1 includeSet2HighpH1 includeSet2PosTHighpH1];
% includeSet2 = [includeSet1HighpH2 includeSet1PosTHighpH2 includeSet2HighpH2];

%%%%%%% High pH  - positive T-wave 

% includeSet1 = [includeSet1HighpH1 includeSet2HighpH1];
% includeSet2 = [includeSet1HighpH2 includeSet2HighpH2];

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


segmentLocation = segmentLocation';
predParam = pH;

featureVectorOriginal = [heartRate' meanTbyQRSRatio' meanTSlope' meanT'];

% featureVectorOriginal = [heartRate' meanTbyQRSRatio' meanTSlope' meanT' ... 
%     normQTOffDuration' normQTDuration' normTTOffDuration' ... 
%     normEnergyPromPCA' normEnergyQTOffPromPCA'];



weights = dec2bin(1:1:2^(length(featureVectorOriginal(1,:)))-1)-'0';

rmsePredParamLinMean = zeros(1,length(weights(:,1)));
rmsePredParamLinMedian = zeros(1,length(weights(:,1)));
rmsePredParamLinVar = zeros(1,length(weights(:,1)));

rmsePredParamSVRMean = zeros(1,length(weights(:,1)));
rmsePredParamSVRMedian = zeros(1,length(weights(:,1)));
rmsePredParamSVRVar = zeros(1,length(weights(:,1)));


predParamEstSVRPlot = zeros(length(weights(:,1)),length(featureVectorOriginal(:,1)));
predParamEstLinPlot = zeros(length(weights(:,1)),length(featureVectorOriginal(:,1)));
predParamTestPlot = zeros(length(weights(:,1)),length(featureVectorOriginal(:,1)));

for featureNum = 1:1:length(weights(:,1))
    featureVector = featureVectorOriginal.*weights(featureNum,:);
    
    featureVector(:,weights(featureNum,:)==0) = [];

    predParamShift = predParam';
    featureVectorShift = featureVector;

    rmsePredParamLin = zeros(1,length(featureVector(:,1)));
    rmsePredParamSVR = zeros(1,length(featureVector(:,1)));
    predictedSegment = zeros(1,length(featureVector(:,1)));

    for segmentNum = 1:1:length(featureVector(:,1))

        predParamTrain = predParamShift(1:end - 1);
        featureVectorTrain = featureVectorShift(1:end - 1,:);

        predParamTest = predParamShift(end);
        featureVectorTest = featureVectorShift(end,:);

        inData = [predParamTrain featureVectorTrain;predParamTest featureVectorTest];
        normArg = 1;
        maxValParamTrain = zeros(1,length(inData(1,:)));
        minValParamTrain = zeros(1,length(inData(1,:)));
        for var2=1:1:length(inData(1,:))
            [inData(:,var2),maxValParamTrain(:,var2),minValParamTrain(:,var2)] = normalizeData(inData(:,var2),normArg);
        end

        predParamTrain = inData(1:end - 1,1);
        featureVectorTrain = inData(1:end - 1,2:end);

        predParamTest = inData(end,1);
        featureVectorTest = inData(end,2:end);


%         svmmodel = svmtrain(predParamTrain,featureVectorTrain,'-s 3 -t 0 -q');  % Low pH - T
%         svmmodel = svmtrain(predParamTrain,featureVectorTrain,'-s 3 -t 3 -q');  % Low pH
%         svmmodel = svmtrain(predParamTrain,featureVectorTrain,'-s 3 -t 2 -q');  % High pH - T 3-3
%         svmmodel = svmtrain(predParamTrain,featureVectorTrain,'-s 4 -t 2 -q');  % High pH 4-1 is also good
%         svmmodel = svmtrain(predParamTrain,featureVectorTrain,'-s 3 -t 2 -q');  % Whole - T
%         svmmodel = svmtrain(predParamTrain,featureVectorTrain,'-s 3 -t 2 -q');  % Whole
        svmmodel = svmtrain(predParamTrain,featureVectorTrain,'-s 3 -t 2 -q');  % Equal - T            
%         svmmodel = svmtrain(predParamTrain,featureVectorTrain,'-s 3 -t 2 -q');  % Equal
        [predClassSVM, ~, ~] = svmpredict(predParamTest,featureVectorTest,svmmodel);
        
        
        predParamEstSVR = predClassSVM;    
        predParamEstSVR = denormalizeData(predParamEstSVR,maxValParamTrain(:,1),minValParamTrain(:,1),normArg);

        coeffsLin = regress(predParamTrain,[ones(length(featureVectorTrain(:,1)),1) featureVectorTrain]);
        predParamEstLin = 0;
        featureVectorTestLin = [ones(length(featureVectorTest(:,1)),1) featureVectorTest];
        for var2=1:1:length(featureVectorTestLin(1,:))
            predParamEstLin = predParamEstLin + coeffsLin(var2)*featureVectorTestLin(:,var2);
        end

        predParamEstLin = predParamEstLin';
        predParamEstLin = denormalizeData(predParamEstLin,maxValParamTrain(:,1),minValParamTrain(:,1),normArg);

        predParamTest = denormalizeData(predParamTest,maxValParamTrain(:,1),minValParamTrain(:,1),normArg);

        predParamEstSVRPlot(featureNum,segmentNum) = predParamEstSVR; 
        predParamEstLinPlot(featureNum,segmentNum) = predParamEstLin;
        predParamTestPlot(featureNum,segmentNum) = predParamTest;


        rmsePredParamLin(segmentNum) = sqrt(sum((predParamTest' - predParamEstLin).^2)/length(predParamTest'));
        rmsePredParamSVR(segmentNum) = sqrt(sum((predParamTest' - predParamEstSVR).^2)/length(predParamTest'));
        predictedSegment(segmentNum) = segmentLocation(end);

        fprintf('RMSE...  Linear:%f, SVR:%f, Predicted Segment :%d \n',rmsePredParamLin(segmentNum),rmsePredParamSVR(segmentNum),predictedSegment(segmentNum));

        featureVectorShift = circshift(featureVectorShift,1,1);
        predParamShift = circshift(predParamShift,1,1);
        segmentLocation = circshift(segmentLocation,1,1);
    end
    
    rmsePredParamLinMean(featureNum) = mean(rmsePredParamLin);
    rmsePredParamLinMedian(featureNum) = median(rmsePredParamLin);
    rmsePredParamLinVar(featureNum) = var(rmsePredParamLin);
    
    rmsePredParamSVRMean(featureNum) = mean(rmsePredParamSVR);
    rmsePredParamSVRMedian(featureNum) = median(rmsePredParamSVR);
    rmsePredParamSVRVar(featureNum) = var(rmsePredParamSVR);
    


    fprintf('\n \n RMSE... \n Linear:%f %f %f \n SVR:%f %f %f \n',mean(rmsePredParamLin),median(rmsePredParamLin),var(rmsePredParamLin), ... 
        mean(rmsePredParamSVR),median(rmsePredParamSVR),var(rmsePredParamSVR));
    fprintf('\n Original:%f %f %f \n Linear:%f %f %f \n SVR:%f %f %f \n',mean(predParamTestPlot),median(predParamTestPlot),var(predParamTestPlot), ... 
        mean(predParamEstLinPlot),median(predParamEstLinPlot),var(predParamEstLinPlot), ... 
        mean(predParamEstSVRPlot),median(predParamEstSVRPlot),var(predParamEstSVRPlot));

end


[rmseValLin,rmseIndLin] = min(rmsePredParamLinMean);
[rmseValSVR,rmseIndSVR] = min(rmsePredParamSVRMean);

barGraphData = [rmsePredParamLinMean' rmsePredParamLinMedian'];
figure
Groups = categorical({'T','TSl','TSl + T','T/QRS','T/QRS + T','T/QRS + TSl','T/QRS + TSl + T','HR', 'HR + T', ... 
    'HR + TSl','HR + TSl + T','HR + T/QRS','HR + T/QRS + T','HR + T/QRS + TSl','HR + T/QRS + TSl + T'});
Groups = reordercats(Groups,{'T','TSl','TSl + T','T/QRS','T/QRS + T','T/QRS + TSl','T/QRS + TSl + T','HR', 'HR + T', ... 
    'HR + TSl','HR + TSl + T','HR + T/QRS','HR + T/QRS + T','HR + T/QRS + TSl','HR + T/QRS + TSl + T'});
bb = bar(Groups,barGraphData);
bb(1).FaceColor = 'r';
bb(2).FaceColor = 'g';
xtickangle(90)
xlabel('Features')
ylabel('RMSE pH')
legend({'Mean', 'Median'},'Location','southeastoutside')



barGraphData = [rmsePredParamSVRMean' rmsePredParamSVRMedian'];
figure
Groups = categorical({'T','TSl','TSl + T','T/QRS','T/QRS + T','T/QRS + TSl','T/QRS + TSl + T','HR', 'HR + T', ... 
    'HR + TSl','HR + TSl + T','HR + T/QRS','HR + T/QRS + T','HR + T/QRS + TSl','HR + T/QRS + TSl + T'});
Groups = reordercats(Groups,{'T','TSl','TSl + T','T/QRS','T/QRS + T','T/QRS + TSl','T/QRS + TSl + T','HR', 'HR + T', ... 
    'HR + TSl','HR + TSl + T','HR + T/QRS','HR + T/QRS + T','HR + T/QRS + TSl','HR + T/QRS + TSl + T'});
bb = bar(Groups,barGraphData);
bb(1).FaceColor = 'r';
bb(2).FaceColor = 'g';
xtickangle(90)
xlabel('Features')
ylabel('RMSE pH')
legend({'Mean', 'Median'},'Location','southeastoutside')


fprintf('\n \n RMSE... \n Linear:%f %f %s \n SVR:%f %f %s \n',rmseValLin,rmsePredParamLinMedian(rmseIndLin),Groups(rmseIndLin),rmseValSVR,rmsePredParamSVRMedian(rmseIndSVR),Groups(rmseIndSVR));
    

figure
stem(predParamTestPlot(rmseIndLin,:),'g', 'filled','MarkerSize',10);
hold on;
stem(predParamEstLinPlot(rmseIndLin,:),'squareb', 'filled','MarkerSize',8);
hold on;
stem(predParamEstSVRPlot(rmseIndSVR,:),':diamondr', 'filled','MarkerSize',6);
hold on;
hold on;
axis([0 length(predParamTestPlot(rmseIndLin,:)) + 1 6.5 8])
xlabel('Segments')
ylabel('pH ')
legend('pH (Original)','estimated pH (LR)','estimated pH (SVR)')




