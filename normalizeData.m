function [normdata,maxVal,minVal] = normalizeData(data,arg)

% arg = 0 [0,1] and arg = 1 for [-1,1] 

maxVal = max(data);
minVal = min(data);

if (maxVal - minVal) == 0
    normdata = (data - minVal);
else
    if arg == 0
        normdata = (data - minVal)./(maxVal - minVal);
    elseif arg == 1
        normdata = -1 + 2.*(data - minVal)./(maxVal - minVal);
    end
end