function [signalBLRemoved] = baselineRemoval(inSignal,sampFreq,order)
% % % First order low pass butterworth filter.

% % % Inputs:
% % % inSignal: Input signal
% % % sampFreq: sampling Frequency in Hz
% % % 
% % % Output:
% % % signalBLRemoved: baseline removed signal

% % % Ref. [1.]  M Varanini, G Tartarisco, L Billeci, A Macerata, G Pioggia
% % % and R Balocchi, "An efficient unsupervised fetal {QRS} complex detection
% % % from abdominal maternal ECG", Physiological Measurement, Vol. 35, No. 8,
% % % 1607-1619, 2014.


if(nargin<3)
    order = 1;  % Assume order 1 if not specified.
end

cutOffFreq = 0.67; %
normCutOffFreq = cutOffFreq/(sampFreq/2); % between [0,1]
[bCoeff,aCoeff]= butter(order,normCutOffFreq,'low');
baseline = filtfilt(bCoeff,aCoeff,inSignal);  %  estimated baseline

%  subtract estimated baseline from original signal
signalBLRemoved = inSignal - baseline;
