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

covTrain = cov(featureVectorTrain);
[eigVecTrain, eigValTrain] = eig(covTrain,'nobalance');
[eigValTrain,orderTrain] = sort(diag(eigValTrain),'descend');  %# sort eigenvalues in descending order
eigVecTrain = eigVecTrain(:,orderTrain);


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

covTest = cov(featureVectorTest);
[eigVecTest, eigValTest] = eig(covTest,'nobalance');
[eigValTest,orderTest] = sort(diag(eigValTest),'descend');  %# sort eigenvalues in descending order
eigVecTest = eigVecTest(:,orderTest);

%% Training

accuracyComp = zeros(1,totalFeatures);
sensitivity = zeros(1,totalFeatures);
specificity = zeros(1,totalFeatures);
geometricMean = zeros(1,totalFeatures);
FScore = zeros(1,totalFeatures);

prominentFeatTrain = featureVectorTrain*eigVecTrain;
prominentFeatTest = featureVectorTest*eigVecTrain;

for selectedFeaturesNo=1:1:totalFeatures

    featureVectorTrain = prominentFeatTrain(:,1:selectedFeaturesNo);
    featureVectorTest = prominentFeatTest(:,1:selectedFeaturesNo);
    
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

    accuracyComp(selectedFeaturesNo) = 100*(truePosCount + trueNegCount)/(truePosCount + trueNegCount + falsePosCount + falseNegCount);
    sensitivity(selectedFeaturesNo) = 100*truePosCount/(truePosCount + falseNegCount);
    specificity(selectedFeaturesNo) = 100*trueNegCount/(trueNegCount + falsePosCount);
    geometricMean(selectedFeaturesNo) = sqrt(sensitivity(selectedFeaturesNo)*specificity(selectedFeaturesNo));
    FScore(selectedFeaturesNo) = 100*(2*truePosCount)/(2*truePosCount + falsePosCount + falseNegCount);
end


figure;
stem(accuracyComp,'b');
hold on;
plot(accuracyComp,'r');
axis([0 25 0 100])
title('PCA-SVC: No. of Features vs Accuracy')
xlabel('Features')
ylabel('Accuracy(%)')
hold on;



