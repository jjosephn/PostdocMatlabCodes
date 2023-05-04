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

includeSet1 = [includeSet1HighpH1 includeSet1PosTHighpH1 includeSet2HighpH1 includeSet2PosTHighpH1 includeSet3LowpH1];
includeSet2 = [includeSet1HighpH2 includeSet1PosTHighpH2 includeSet2HighpH2 includeSet1LowpH2 includeSet1PosTLowpH2 includeSet2LowpH2 includeSet3LowpH2 includeSet3PosTLowpH2];

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

% includeSet1 = [145 146 ... 
%     2 35 43 94 230];
% includeSet2 = [162 166 173 174 175 ... 
%     27 40];

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
% featureVector = [ones(length(meanTbyQRSRatio'),1) meanTbyQRSRatio'];
% featureVector = [ones(length(RRAvg'),1) RRAvg'];
% featureVector = [ones(length(meanT'),1) meanT'];
featureVector = [meanTbyQRSRatio' RRAvg' meanT' meanTSlope' heartRate'];


g1 = find(pH<=7.2);
g2 = find(predParam>7.2&predParam<7.25);
g3 = find(pH>=7.25);

%%

normArg = 1;
maxValParamTrain = zeros(1,length(featureVector(1,:)));
minValParamTrain = zeros(1,length(featureVector(1,:)));
for var2=1:1:length(featureVector(1,:))
    [featureVector(:,var2),maxValParamTrain(:,var2),minValParamTrain(:,var2)] = normalizeData(featureVector(:,var2),normArg);
end

meanTbyQRSRatio = featureVector(:,1)';
RRAvg = featureVector(:,2)';
meanT = featureVector(:,3)';
meanTSlope = featureVector(:,4)';
heartRate = featureVector(:,5)';

xx1g1 = [ones(length(meanTbyQRSRatio(g1)'),1) meanTbyQRSRatio(g1)'];
xx2g1 = [ones(length(RRAvg(g1)'),1) RRAvg(g1)'];
xx3g1 = [ones(length(meanT(g1)'),1) meanT(g1)'];
xx4g1 = [ones(length(meanTSlope(g1)'),1) meanTSlope(g1)'];
xx5g1 = [ones(length(heartRate(g1)'),1) heartRate(g1)'];

coeffs1g1 =regress(pH(g1)',xx1g1);
coeffs2g1 =regress(pH(g1)',xx2g1);
coeffs3g1 =regress(pH(g1)',xx3g1);
coeffs4g1 =regress(pH(g1)',xx4g1);
coeffs5g1 =regress(pH(g1)',xx5g1);
pHest1g1 = coeffs1g1(1) + coeffs1g1(2).*meanTbyQRSRatio(g1)';
pHest2g1 = coeffs2g1(1) + coeffs2g1(2).*RRAvg(g1)';
pHest3g1 = coeffs3g1(1) + coeffs3g1(2).*meanT(g1)';
pHest4g1 = coeffs4g1(1) + coeffs4g1(2).*meanTSlope(g1)';
pHest5g1 = coeffs5g1(1) + coeffs5g1(2).*heartRate(g1)';


cc1g1 = corrcoef(meanTbyQRSRatio(g1),pH(g1));
cc2g1 = corrcoef(RRAvg(g1),pH(g1));
cc3g1 = corrcoef(meanT(g1),pH(g1));
cc4g1 = corrcoef(meanTSlope(g1),pH(g1));
cc5g1 = corrcoef(heartRate(g1),pH(g1));

fprintf('T/QRS: %f... RRAvg: %f ... T: %f ... TSlope: %f ... HR: %f \n',cc1g1(1,2),cc2g1(1,2),cc3g1(1,2),cc4g1(1,2),cc5g1(1,2));


xx1g2 = [ones(length(meanTbyQRSRatio(g2)'),1) meanTbyQRSRatio(g2)'];
xx2g2 = [ones(length(RRAvg(g2)'),1) RRAvg(g2)'];
xx3g2 = [ones(length(meanT(g2)'),1) meanT(g2)'];
xx4g2 = [ones(length(meanTSlope(g2)'),1) meanTSlope(g2)'];
xx5g2 = [ones(length(heartRate(g2)'),1) heartRate(g2)'];

coeffs1g2 =regress(pH(g2)',xx1g2);
coeffs2g2 =regress(pH(g2)',xx2g2);
coeffs3g2 =regress(pH(g2)',xx3g2);
coeffs4g2 =regress(pH(g2)',xx4g2);
coeffs5g2 =regress(pH(g2)',xx5g2);
pHest1g2 = coeffs1g2(1) + coeffs1g2(2).*meanTbyQRSRatio(g2)';
pHest2g2 = coeffs2g2(1) + coeffs2g2(2).*RRAvg(g2)';
pHest3g2 = coeffs3g2(1) + coeffs3g2(2).*meanT(g2)';
pHest4g2 = coeffs4g2(1) + coeffs4g2(2).*meanTSlope(g2)';
pHest5g2 = coeffs5g2(1) + coeffs5g2(2).*heartRate(g2)';


cc1g2 = corrcoef(meanTbyQRSRatio(g2),pH(g2));
cc2g2 = corrcoef(RRAvg(g2),pH(g2));
cc3g2 = corrcoef(meanT(g2),pH(g2));
cc4g2 = corrcoef(meanTSlope(g2),pH(g2));
cc5g2 = corrcoef(heartRate(g2),pH(g2));

fprintf('T/QRS: %f... RRAvg: %f ... T: %f ... TSlope: %f ... HR: %f \n',cc1g2(1,2),cc2g2(1,2),cc3g2(1,2),cc4g2(1,2),cc5g2(1,2));



xx1g3 = [ones(length(meanTbyQRSRatio(g3)'),1) meanTbyQRSRatio(g3)'];
xx2g3 = [ones(length(RRAvg(g3)'),1) RRAvg(g3)'];
xx3g3 = [ones(length(meanT(g3)'),1) meanT(g3)'];
xx4g3 = [ones(length(meanTSlope(g3)'),1) meanTSlope(g3)'];
xx5g3 = [ones(length(heartRate(g3)'),1) heartRate(g3)'];

coeffs1g3 =regress(pH(g3)',xx1g3);
coeffs2g3 =regress(pH(g3)',xx2g3);
coeffs3g3 =regress(pH(g3)',xx3g3);
coeffs4g3 =regress(pH(g3)',xx4g3);
coeffs5g3 =regress(pH(g3)',xx5g3);
pHest1g3 = coeffs1g3(1) + coeffs1g3(2).*meanTbyQRSRatio(g3)';
pHest2g3 = coeffs2g3(1) + coeffs2g3(2).*RRAvg(g3)';
pHest3g3 = coeffs3g3(1) + coeffs3g3(2).*meanT(g3)';
pHest4g3 = coeffs4g3(1) + coeffs4g3(2).*meanTSlope(g3)';
pHest5g3 = coeffs5g3(1) + coeffs5g3(2).*heartRate(g3)';


cc1g3 = corrcoef(meanTbyQRSRatio(g3),pH(g3));
cc2g3 = corrcoef(RRAvg(g3),pH(g3));
cc3g3 = corrcoef(meanT(g3),pH(g3));
cc4g3 = corrcoef(meanTSlope(g3),pH(g3));
cc5g3 = corrcoef(heartRate(g3),pH(g3));

fprintf('T/QRS: %f... RRAvg: %f ... T: %f ... TSlope: %f ... HR: %f \n',cc1g3(1,2),cc2g3(1,2),cc3g3(1,2),cc4g3(1,2),cc5g3(1,2));

%% Plots

group = zeros(1,length(pH));
group(pH<=7.2) = 1;
group(pH>7.2&predParam<7.25) = 2;
group(pH>=7.25) = 3;


figure;
gscatter(RRAvg,pH,group,'rgb','xo*')
hold on;
plot(RRAvg(g1),pHest2g1,'r')
hold on;
plot(RRAvg(g2),pHest2g2,'g')
hold on;
plot(RRAvg(g3),pHest2g3,'b')
% title('RR Average vs pH')
% axis([0 length(predParamIndex) + 1 6.85 7.85])
xlabel('RR Average')
ylabel('pH ')
% legend('pH > 7.2','7.2< pH<7.25','pH >= 7.25','fitted line:pH > 7.2','fitted line:7.2< pH<7.25 ','fitted line: pH >= 7.25')
% legend('boxoff')
legend('Off')


figure;
gscatter(meanTbyQRSRatio,pH,group,'rgb','xo*')
hold on;
plot(meanTbyQRSRatio(g1),pHest1g1,'r')
hold on;
plot(meanTbyQRSRatio(g2),pHest1g2,'g')
hold on;
plot(meanTbyQRSRatio(g3),pHest1g3,'b')
% title('T/QRS vs pH')
% axis([0 length(predParamIndex) + 1 6.85 7.85])
xlabel('T/QRS Ratio')
ylabel('pH ')
legend('pH > 7.2','7.2< pH<7.25','pH >= 7.25','fitted line:pH > 7.2','fitted line:7.2< pH<7.25 ','fitted line: pH >= 7.25')
legend('boxoff')

figure;
gscatter(meanT,pH,group,'rgb','xo*')
hold on;
plot(meanT(g1),pHest3g1,'r')
hold on;
plot(meanT(g2),pHest3g2,'g')
hold on;
plot(meanT(g3),pHest3g3,'b')
% title('T vs pH')
% axis([0 length(predParamIndex) + 1 6.85 7.85])
xlabel('T Amplitude')
ylabel('pH ')
% legend('pH > 7.2','7.2< pH<7.25','pH >= 7.25','fitted line:pH > 7.2','fitted line:7.2< pH<7.25 ','fitted line: pH >= 7.25')
% legend('boxoff')
legend('Off')

figure;
gscatter(meanTSlope,pH,group,'rgb','xo*')
hold on;
plot(meanTSlope(g1),pHest4g1,'r')
hold on;
plot(meanTSlope(g2),pHest4g2,'g')
hold on;
plot(meanTSlope(g3),pHest4g3,'b')
% title('T-TOff Slope vs pH')
% axis([0 length(predParamIndex) + 1 6.85 7.85])
xlabel('T Slope')
ylabel('pH ')
% legend('pH > 7.2','7.2< pH<7.25','pH >= 7.25','fitted line:pH > 7.2','fitted line:7.2< pH<7.25 ','fitted line: pH >= 7.25')
% legend('boxoff')
legend('Off')

figure;
gscatter(heartRate,pH,group,'rgb','xo*')
hold on;
plot(heartRate(g1),pHest5g1,'r')
hold on;
plot(heartRate(g2),pHest5g2,'g')
hold on;
plot(heartRate(g3),pHest5g3,'b')
% title('Heart Rate vs pH')
% axis([0 length(predParamIndex) + 1 6.85 7.85])
xlabel('Heart Rate')
ylabel('pH ')
% legend('pH > 7.2','7.2< pH<7.25','pH >= 7.25','fitted line:pH > 7.2','fitted line:7.2< pH<7.25 ','fitted line: pH >= 7.25')
% legend('boxoff')
legend('Off')


