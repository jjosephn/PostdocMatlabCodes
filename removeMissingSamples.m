function [rmMissing] = removeMissingSamples(tempFullSignalStage1,tempFullSignalStage2,selectFlag,sampFreq,secondStageStart)

% 0 removes missing samples more than 20 seconds in FHR only. 1 removes missing samples more than 20 seconds from both FHR and UC.

%% Stage I FHR

segLenSeconds = 20;
skipSegmentSamples = segLenSeconds*sampFreq;

tempFullSignalStage11 = NaN*ones(size(tempFullSignalStage1));
var1 = 1;
var2 = 1;
var3 = 1;
[missingSamplesStartFHR1,missingSamplesStopFHR1] = missingSamples(tempFullSignalStage1(1,:));
while var1 <= length(tempFullSignalStage1(1,:))
    if var2 <= length(missingSamplesStartFHR1) && missingSamplesStartFHR1(var2) == var1
        if (missingSamplesStopFHR1(var2) - missingSamplesStartFHR1(var2)) >= skipSegmentSamples
            var1 = missingSamplesStopFHR1(var2) + 1;
            var2 = var2 + 1;
        else
            tempFullSignalStage11(:,var3) = tempFullSignalStage1(:,var1);
            var3 = var3 + 1;
            var2 = var2 + 1;
        end
    else
        tempFullSignalStage11(:,var3) = tempFullSignalStage1(:,var1);
        var3 = var3 + 1;
    end
    var1 = var1 + 1;
end

tempFullSignalStage11(:,isnan(tempFullSignalStage11(1,:))) = [];
rmMissing.tempFullSignalStage11 = tempFullSignalStage11;

%% Stage I UC


tempFullSignalStage12 = NaN*ones(size(tempFullSignalStage11));
var1 = 1;
var2 = 1;
var3 = 1;
[missingSamplesStartUC1,missingSamplesStopUC1] = missingSamples(tempFullSignalStage11(2,:));
while var1 <= length(tempFullSignalStage11(1,:))
    if var2 <= length(missingSamplesStartUC1) && missingSamplesStartUC1(var2) == var1
        if (missingSamplesStopUC1(var2) - missingSamplesStartUC1(var2)) >= skipSegmentSamples
            var1 = missingSamplesStopUC1(var2) + 1;
            var2 = var2 + 1;
        else
            tempFullSignalStage12(:,var3) = tempFullSignalStage11(:,var1);
            var3 = var3 + 1;
            var2 = var2 + 1;
        end
    else
        tempFullSignalStage12(:,var3) = tempFullSignalStage11(:,var1);
        var3 = var3 + 1;
    end
    var1 = var1 + 1;
end
if selectFlag == 0
    tempFullSignalStage12 = tempFullSignalStage11;
else
    tempFullSignalStage12(:,isnan(tempFullSignalStage12(2,:))) = [];
end
rmMissing.tempFullSignalStage12 = tempFullSignalStage12;

if secondStageStart > -1

%% Stage II FHR

    tempFullSignalStage21 = NaN*ones(size(tempFullSignalStage2));
    var1 = 1;
    var2 = 1;
    var3 = 1;
    [missingSamplesStartFHR2,missingSamplesStopFHR2] = missingSamples(tempFullSignalStage2(1,:));
    while var1 <= length(tempFullSignalStage2(1,:))
        if var2 <= length(missingSamplesStartFHR2) && missingSamplesStartFHR2(var2) == var1
            if (missingSamplesStopFHR2(var2) - missingSamplesStartFHR2(var2)) >= skipSegmentSamples
                var1 = missingSamplesStopFHR2(var2) + 1;
                var2 = var2 + 1;
            else
                tempFullSignalStage21(:,var3) = tempFullSignalStage2(:,var1);
                var3 = var3 + 1;
                var2 = var2 + 1;
            end
        else
            tempFullSignalStage21(:,var3) = tempFullSignalStage2(:,var1);
            var3 = var3 + 1;
        end
        var1 = var1 + 1;
    end
    tempFullSignalStage21(:,isnan(tempFullSignalStage21(1,:))) = [];
    rmMissing.tempFullSignalStage21 = tempFullSignalStage21;

%% Stage II UC

    tempFullSignalStage22 = NaN*ones(size(tempFullSignalStage21));
    var1 = 1;
    var2 = 1;
    var3 = 1;
    [missingSamplesStartUC2,missingSamplesStopUC2] = missingSamples(tempFullSignalStage21(2,:));
    while var1 <= length(tempFullSignalStage21(1,:))
        if var2 <= length(missingSamplesStartUC2) && missingSamplesStartUC2(var2) == var1
            if (missingSamplesStopUC2(var2) - missingSamplesStartUC2(var2)) >= skipSegmentSamples
                var1 = missingSamplesStopUC2(var2) + 1;
                var2 = var2 + 1;
            else
                tempFullSignalStage22(:,var3) = tempFullSignalStage21(:,var1);
                var3 = var3 + 1;
                var2 = var2 + 1;
            end
        else
            tempFullSignalStage22(:,var3) = tempFullSignalStage21(:,var1);
            var3 = var3 + 1;
        end
        var1 = var1 + 1;
    end
    if selectFlag == 0
        tempFullSignalStage22 = tempFullSignalStage21;
    else
        tempFullSignalStage22(:,isnan(tempFullSignalStage22(2,:))) = [];
    end
    rmMissing.tempFullSignalStage22 = tempFullSignalStage22;
else
    tempFullSignalStage21 = zeros(2,1);
    tempFullSignalStage22 = zeros(2,1);
    rmMissing.tempFullSignalStage21 = tempFullSignalStage21;
    rmMissing.tempFullSignalStage22 = tempFullSignalStage22;
end