function [UCSignalInfo] = UCInfoExtract(UCSignal,sampFreq)

UCSignal = UCSignal/abs(max(UCSignal));

%% Band Pass filter

cutOffFreqHigh = 1/60;                                         % minimum frequency of occurence =  1 minute (1/60 )
butterOrder = 3;
normCutOffFreq = cutOffFreqHigh*2/sampFreq;                % between [0,1]
[bCoeff,aCoeff]= butter(butterOrder,normCutOffFreq);
UCBandPassSignal = filtfilt(bCoeff,aCoeff,UCSignal);                         % estimated baseline
UCBandPassSignal = UCBandPassSignal/max(abs(UCBandPassSignal));

%% Find Peaks

UCSignalMinDistSeconds = 40;
UCSignalMinDistSamples = UCSignalMinDistSeconds*sampFreq;
conractionFlag = 0;
if length(UCBandPassSignal) > 2*UCSignalMinDistSamples
    [UCBandPassSignalPeaks,UCBandPassSignalLocations,~,~] = findpeaks(UCBandPassSignal,'MinPeakDistance',UCSignalMinDistSamples);

    UCBandPassSignalLocations(UCBandPassSignalPeaks<0.4*mean(UCBandPassSignalPeaks)) = [];

    valleyLocs = zeros(1,length(UCBandPassSignalLocations));
    valleyVal = zeros(1,length(UCBandPassSignalLocations));
    for peakNos = 1:1:length(UCBandPassSignalLocations)
        if peakNos == 1
            [valleyVal(peakNos),valleyLocs(peakNos)] = min(UCBandPassSignal(1:UCBandPassSignalLocations(peakNos)));
        else
            [valleyVal(peakNos),valleyLocs(peakNos)] = min(UCBandPassSignal(UCBandPassSignalLocations(peakNos - 1):UCBandPassSignalLocations(peakNos)));
            valleyLocs(peakNos) = UCBandPassSignalLocations(peakNos - 1) + valleyLocs(peakNos) - 1;
            if peakNos == length(UCBandPassSignalLocations)
                [valleyVal(peakNos + 1),valleyLocs(peakNos + 1)] = min(UCBandPassSignal(UCBandPassSignalLocations(peakNos):end));
                valleyLocs(peakNos + 1) = UCBandPassSignalLocations(peakNos) + valleyLocs(peakNos + 1) - 1;
            end
        end
    end

    UCBandPassStartLocation = zeros(1,length(UCBandPassSignalLocations));
    UCBandPassStopLocation = zeros(1,length(UCBandPassSignalLocations));
    UCBandPassStartVal = zeros(1,length(UCBandPassSignalLocations));
    UCBandPassStopVal = zeros(1,length(UCBandPassSignalLocations));
    for peakVar1 = 1:1:length(UCBandPassSignalLocations)
        if (length(UCBandPassSignalLocations) < length(valleyVal))
            UCBandPassStartLocation(peakVar1) = valleyLocs(peakVar1);
            UCBandPassStartVal(peakVar1) = valleyVal(peakVar1);
            UCBandPassStopLocation(peakVar1) = valleyLocs(peakVar1 + 1);
            UCBandPassStopVal(peakVar1) = valleyVal(peakVar1 + 1);
            conractionFlag = 1;
        end
    end
    UCSignalStartLocation = zeros(1,length(UCBandPassSignalLocations));
    UCSignalStopLocation = zeros(1,length(UCBandPassSignalLocations));
    UCSignalStartVal = zeros(1,length(UCBandPassSignalLocations));
    UCSignalStopVal = zeros(1,length(UCBandPassSignalLocations));

end

if conractionFlag == 1
    searchPercentage = 0.6;
    slopeWindowSizePercentage = 0.2;
    for peakVar1=1:1:length(UCBandPassStartLocation)
        searchOnStart = UCBandPassStartLocation(peakVar1);
        searchOnStop = round((UCBandPassSignalLocations(peakVar1) - UCBandPassStartLocation(peakVar1))*searchPercentage);
        searchOnStop = searchOnStop + searchOnStart - 1;
        if searchOnStop >= length(UCSignal)
            searchOnStop = length(UCSignal);
        end
        peakVar3 = 1;
        slopeWindowSize = round((UCBandPassSignalLocations(peakVar1) - UCBandPassStartLocation(peakVar1))*slopeWindowSizePercentage);
        startSlope = zeros(1,searchOnStop - searchOnStart - slopeWindowSize + 2);
        for peakVar2 = searchOnStart:1:searchOnStop - slopeWindowSize + 1
            startSlope(peakVar3) = mean(diff(UCSignal(peakVar2:peakVar2 + slopeWindowSize - 1)));
            peakVar3 = peakVar3 + 1;
        end
        [~, tempIndx1] = max(startSlope);
        UCSignalStartLocation(peakVar1) = searchOnStart + tempIndx1;
        UCSignalStartVal(peakVar1) = UCSignal(UCSignalStartLocation(peakVar1));

        searchOffStart = round((UCBandPassStopLocation(peakVar1) - UCBandPassSignalLocations(peakVar1))*searchPercentage);                
        searchOffStart = UCBandPassStopLocation(peakVar1) - searchOffStart + 1;
        searchOffStop = UCBandPassStopLocation(peakVar1);

        peakVar3 = 1;
        slopeWindowSize = round((UCBandPassStopLocation(peakVar1) - UCBandPassSignalLocations(peakVar1))*slopeWindowSizePercentage);
        stopSlope = zeros(1,searchOffStop - searchOffStart - slopeWindowSize + 2);   
        flipSignal = fliplr(UCSignal(searchOffStart:searchOffStop));
        for peakVar2 = 1:1:length(flipSignal) - slopeWindowSize + 1
            stopSlope(peakVar3) = mean(diff(flipSignal(peakVar2:peakVar2 + slopeWindowSize - 1)));
            peakVar3 = peakVar3 + 1;
        end
        [~, tempIndx2] = max(stopSlope);
        UCSignalStopLocation(peakVar1) = searchOffStop - tempIndx2;
        UCSignalStopVal(peakVar1) = UCSignal(UCSignalStopLocation(peakVar1));
    end

    UCSignalPeaks = zeros(1,length(UCBandPassSignalLocations));
    UCSignalLocations = zeros(1,length(UCBandPassSignalLocations));            
    for peakNos = 1:1:length(UCBandPassSignalLocations)
        [UCSignalPeaks(peakNos),UCSignalLocations(peakNos)] = max(UCSignal(UCSignalStartLocation(peakNos):UCSignalStopLocation(peakNos)));
        UCSignalLocations(peakNos) = UCSignalStartLocation(peakNos) + UCSignalLocations(peakNos) - 1;
    end
    UCSignalInfo.UCSignalLocations = UCSignalLocations;
    UCSignalInfo.UCSignalPeaks = UCSignalPeaks;
    UCSignalInfo.UCSignalPeaksNumbers = length(UCSignalPeaks);
    UCSignalInfo.UCSignalStartLocation = UCSignalStartLocation;
    UCSignalInfo.UCSignalStartVal = UCSignalStartVal;
    UCSignalInfo.UCSignalStopLocation = UCSignalStopLocation;
    UCSignalInfo.UCSignalStopVal = UCSignalStopVal;
else
    UCSignalInfo.UCSignalLocations = [];
    UCSignalInfo.UCSignalPeaks = [];
    UCSignalInfo.UCSignalPeaksNumbers = 0;
    UCSignalInfo.UCSignalStartLocation = [];
    UCSignalInfo.UCSignalStartVal = [];
    UCSignalInfo.UCSignalStopLocation = [];
    UCSignalInfo.UCSignalStopVal = [];
end
