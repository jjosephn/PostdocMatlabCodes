function [signalNoiseRemoved] = highFreqNoiseRemoval(inSignal,sampFreq,butterOrder,cutOffFreqLow)
% % % First order low pass butterworth filter.

% % % Inputs:
% % % inSignal: Input signal
% % % sampFreq: sampling Frequency in Hz
% % % butterOrder: order
% % % cutOffFreqLow: cut-off frequency
% % % 
% % % Output:
% % % signalBLRemoved: baseline removed signal

if(nargin<3)
    butterOrder = 1;  % Assume order 1 if not specified.
end
if(nargin<4)
    cutOffFreqLow = 50; % assume 50 Hz if not specified
end

%% Low pass filter

% cutOffFreqLow = 50; %
% butterOrder = 1;
normCutOffFreq = cutOffFreqLow/(sampFreq/2); % between [0,1]
[bCoeff,aCoeff]= butter(butterOrder,normCutOffFreq,'low');
signalLowPass = filtfilt(bCoeff,aCoeff,inSignal);                         % estimated baseline
signalNoiseRemoved = signalLowPass;


