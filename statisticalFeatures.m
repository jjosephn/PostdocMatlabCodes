function [statFeatures] = statisticalFeatures(signal)

statFeatures.statFeatMean = mean(signal);
statFeatures.statFeatMedian = median(signal);
statFeatures.statFeatSD = std(signal);
statFeatures.statFeatSE = statFeatures.statFeatSD/sqrt(length(signal));
statFeatures.statFeatVar = var(signal);
statFeatures.statFeatMax = max(signal);
statFeatures.statFeatMin = min(signal);
statFeatures.statFeatMaxMinDiff = statFeatures.statFeatMax - statFeatures.statFeatMin;
statFeatures.statFeatMeanAbsDev = mad(signal,0);
statFeatures.statFeatMedianAbsDev = mad(signal,1);