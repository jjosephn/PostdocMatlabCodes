function [missingSamplesStart,missingSamplesStop] = missingSamples(signal)

missingSamplesStart = zeros(1,500);
missingSamplesStop = zeros(1,500);
signalLength = length(signal);
sample = 1;
var1 = 1;
while (sample <= signalLength)
    if signal(sample) == 0
        missingSamplesStart(var1) = sample;
        flag = 1;
        while flag
            if signal(sample) > 0
                missingSamplesStop(var1) = sample - 1;
                flag = 0;
            elseif sample == signalLength
                missingSamplesStop(var1) = sample;
                flag = 0;
            end
            sample = sample + 1;
        end
        var1 = var1 + 1;
    else
        sample = sample + 1;
    end
end

missingSamplesStart(missingSamplesStart == 0) = [];
missingSamplesStop(missingSamplesStop == 0) = [];

