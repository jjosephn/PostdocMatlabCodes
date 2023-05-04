function [orgdata] = denormalizeData(normdata,maxVal,minVal,arg)

% arg = 0 [0,1] and arg = 1 for [-1,1] 

if arg == 0
    orgdata = (normdata.*(maxVal - minVal)) + minVal;
elseif arg == 1
    orgdata = (normdata.*(maxVal - minVal) + maxVal + minVal)./2;
end