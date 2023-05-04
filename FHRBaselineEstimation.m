function [FHRBaselineSegment, continuousSegments] = FHRBaselineEstimation(FHRSignal)

virtualBaseline = mean(FHRSignal);
continuousSegments = NaN*ones(1,length(FHRSignal));
alpha = 8;

for var1 = 1:1:length(FHRSignal)
    if (FHRSignal(var1) <= virtualBaseline + alpha) && (FHRSignal(var1) >= virtualBaseline - alpha)
        continuousSegments(var1) = FHRSignal(var1);
    else
        continuousSegments(var1) = NaN;
    end        
end

continuousSegments(:,isnan(continuousSegments)) = [];
FHRBaselineSegment = mean(continuousSegments);
