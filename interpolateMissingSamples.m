function [interMissing] = interpolateMissingSamples(rmMissing,secondStageStart)

tempFullSignalStage12 = rmMissing.tempFullSignalStage12;
tempFullSignalStage22 = rmMissing.tempFullSignalStage22;


interType = 'spline';    

InterFHR1 = tempFullSignalStage12(1,:);
InterFHR1(InterFHR1==0) = NaN;
InterFHR1 = fillmissing(InterFHR1,interType,'SamplePoints',1:length(InterFHR1));

InterUC1 = tempFullSignalStage12(2,:);
InterUC1(InterUC1==0) = NaN;
InterUC1 = fillmissing(InterUC1,interType,'SamplePoints',1:length(InterUC1));

interMissing.fullSignalStageI = [InterFHR1;InterUC1];


if secondStageStart > -1        
    InterFHR2 = tempFullSignalStage22(1,:);
    InterFHR2(InterFHR2==0) = NaN;
    InterFHR2 = fillmissing(InterFHR2,interType,'SamplePoints',1:length(InterFHR2));

    InterUC2 = tempFullSignalStage22(2,:);
    InterUC2(InterUC2==0) = NaN;
    InterUC2 = fillmissing(InterUC2,interType,'SamplePoints',1:length(InterUC2));

    interMissing.fullSignalStageII = [InterFHR2;InterUC2];

else
    interMissing.fullSignalStageII = zeros(2,1);
end