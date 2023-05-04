clc;
close all;
clearvars;

[parameters,signal,info] = initDatabase();

pHCutOff = 7.20;

load 'featuresExtractionV1o010MtsSegmentEqual';

segmentNo = 1;
featStart = 1;
featStop = length(features.train.(features.train.fileInfo{1}).featureVectorTrain);
normArg = 1;
totalFeatures = featStop - featStart + 1;
totalParams = length(features.train.(features.train.fileInfo{1}).clinicalParamsVectorTrain);

clinicalParamsVectorTrain = zeros(length(features.train.fileInfo),totalParams);
featureVectorTrain = zeros(length(features.train.fileInfo),totalFeatures);

clinicalParamsVectorTest = zeros(length(features.test.fileInfo),totalParams);
featureVectorTest = zeros(length(features.test.fileInfo),totalFeatures);

for var1 = 1:1:length(features.train.fileInfo)
    fileValsChar = features.train.fileInfo{var1};
    clinicalParamsVectorTrain(var1,:) = features.train.(fileValsChar).clinicalParamsVectorTrain;
    featureVectorTrain(var1,1:totalFeatures) = features.train.(fileValsChar).featureVectorTrain(featStart:featStop,segmentNo)';
end

for var1 = 1:1:length(features.test.fileInfo)
    fileValsChar = features.test.fileInfo{var1};
    clinicalParamsVectorTest(var1,:) = features.test.(fileValsChar).clinicalParamsVectorTest;
    featureVectorTest(var1,1:totalFeatures) = features.test.(fileValsChar).featureVectorTest(featStart:featStop,segmentNo)';
end

%% Train Data

index = zeros(1,8);
var2 = 1;
for var1 = 1:1:length(featureVectorTrain(:,1))
    if ~isempty(find(isnan(featureVectorTrain(var1,:)), 1))
        index(var2) = var1;
        var2 = var2 + 1;
    elseif ~isempty(find(isnan(clinicalParamsVectorTrain(var1,:)), 1))
        index(var2) = var1;
        var2 = var2 + 1;
    end
end

index(index==0) = [];

featureVectorTrain(index,:) = [];
clinicalParamsVectorTrain(index,:) = [];

clinicalParamsVectorTrainLabel = zeros(length(clinicalParamsVectorTrain(:,1)),1);
for var1 = 1:1:length(clinicalParamsVectorTrain(:,1))
    if clinicalParamsVectorTrain(var1,1) <= pHCutOff
        clinicalParamsVectorTrainLabel(var1,1) = 1;
    else
        clinicalParamsVectorTrainLabel(var1,1) = -1;
    end    
end

for var1=1:1:length(featureVectorTrain(1,:))
    [featureVectorTrain(:,var1),~,~] = normalizeData(featureVectorTrain(:,var1),normArg);
end

%% Test Data

index = zeros(1,3);
var2 = 1;
for var1 = 1:1:length(featureVectorTest(:,1))
    if ~isempty(find(isnan(featureVectorTest(var1,:)), 1))
        index(var2) = var1;
        var2 = var2 + 1;
    elseif ~isempty(find(isnan(clinicalParamsVectorTest(var1,:)), 1))
        index(var2) = var1;
        var2 = var2 + 1;
    end
end

index(index==0) = [];

featureVectorTest(index,:) = [];
clinicalParamsVectorTest(index,:) = [];

for var1=1:1:length(featureVectorTest(1,:))
    [featureVectorTest(:,var1),~,~] = normalizeData(featureVectorTest(:,var1),normArg);
end

clinicalParamsVectorTestLabel = zeros(length(clinicalParamsVectorTest(:,1)),1);
for var1 = 1:1:length(clinicalParamsVectorTest(:,1))
    if clinicalParamsVectorTest(var1,1) <= pHCutOff
        clinicalParamsVectorTestLabel(var1,1) = 1;
    else
        clinicalParamsVectorTestLabel(var1,1) = -1;
    end    
end

%% Training

for var1=1:1:length(featureVectorTrain(1,:))
    [featureVectorTrain(:,var1),~,~] = normalizeData(featureVectorTrain(:,var1),1);
end

for var1=1:1:length(featureVectorTest(1,:))
    [featureVectorTest(:,var1),~,~] = normalizeData(featureVectorTest(:,var1),1);
end

randIndex = randperm(size(featureVectorTrain,1));

featureVectorTrain = featureVectorTrain(randIndex,:);
clinicalParamsVectorTrainLabel = clinicalParamsVectorTrainLabel(randIndex,:);

model = svmtrain(clinicalParamsVectorTrainLabel,featureVectorTrain,'-s 0 -t 0');
[predClinicalParamsVectorLabel, accuracy, probEstimates] = svmpredict(clinicalParamsVectorTestLabel,featureVectorTest,model);

truePosCount = 0;
trueNegCount = 0;
falsePosCount = 0;
falseNegCount = 0;
for var1=1:1:length(predClinicalParamsVectorLabel)
    if clinicalParamsVectorTestLabel(var1) == 1 && predClinicalParamsVectorLabel(var1) == 1
        truePosCount = truePosCount + 1;
    elseif clinicalParamsVectorTestLabel(var1) == -1 && predClinicalParamsVectorLabel(var1) == -1
        trueNegCount = trueNegCount + 1;
    elseif clinicalParamsVectorTestLabel(var1) == -1 && predClinicalParamsVectorLabel(var1) == 1
        falsePosCount = falsePosCount + 1;
    elseif clinicalParamsVectorTestLabel(var1) == 1 && predClinicalParamsVectorLabel(var1) == -1
        falseNegCount = falseNegCount + 1;
    end
end

accuracy = 100*(truePosCount + trueNegCount)/(truePosCount + trueNegCount + falsePosCount + falseNegCount);
sensitivity = 100*truePosCount/(truePosCount + falseNegCount);
specificity = 100*trueNegCount/(trueNegCount + falsePosCount);
geometricMean = sqrt(sensitivity*specificity);
FScore = 100*(2*truePosCount)/(2*truePosCount + falsePosCount + falseNegCount);

fprintf('accuracy = %2.2f %% \n sensitivity = %2.2f %% \n specificity = %2.2f %% \n geometricMean = %2.2f %% \n FScore = %2.2f %% \n', ... 
    accuracy,sensitivity,specificity,geometricMean,FScore);


