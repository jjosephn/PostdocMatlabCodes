function [QRSFeatures]=QRSFeaturesNeo(inSignal,sampFreq)

signalRFeatures = panTompkinNeo(inSignal,sampFreq);

% signalBandPass = signalRFeatures.signal;
% peaksSignalBandPass = signalRFeatures.peaksValueSignal;
% peaksLocationBandPass = signalRFeatures.peaksLocationSignal;

QRSFeatures.signalBandPass = signalRFeatures.signal;
QRSFeatures.peaksSignalBandPass = signalRFeatures.peaksValueSignal;
QRSFeatures.peaksLocationBandPass = signalRFeatures.peaksLocationSignal;
QRSFeatures.RRAverageBandPass1 = signalRFeatures.RRAverageBandPass1;


% 
% %% 1 -40 Hz filter
% 
% inSignal = signalBandPass;
% 
% inSignal = inSignal - mean(inSignal);
% inSignal = inSignal/max(abs(inSignal));
% 
% QRSInitialWidth = 0.080;                                                   % Initial QRS width in ms for fetal ECG (80 -120 ms ???).
% windowLengthQRS = round(QRSInitialWidth*sampFreq);
% if mod(windowLengthQRS,2) == 1
%     windowLengthQRS = windowLengthQRS + 1;
% end
% 
% paddedWindowLength = windowLengthQRS;
% paddedInSignal = [zeros(1,paddedWindowLength) inSignal];
% 
% butterOrder40 = 4;
% cutOffFreqLow40 = 1;
% cutOffFreqHigh40 = 40;
% normCutOffFreq40 = [cutOffFreqLow40 cutOffFreqHigh40]*2/sampFreq;                % between [0,1]
% [bCoeff40,aCoeff40]= butter(butterOrder40,normCutOffFreq40);
% bandPassSignal40 = filtfilt(bCoeff40,aCoeff40,paddedInSignal);
% bandPassSignal40 = bandPassSignal40(paddedWindowLength + 1:end);
% 
% bandPassSignal40 = bandPassSignal40 - mean(bandPassSignal40);
% bandPassSignal40 = bandPassSignal40/max(abs(bandPassSignal40));
% 
% derivativeBPSignal40 = diff(bandPassSignal40)/2;
% 
% derivativeBPSignal40 = derivativeBPSignal40 - mean(derivativeBPSignal40);
% derivativeBPSignal40 = derivativeBPSignal40/max(abs(derivativeBPSignal40));
% 
% derivativeBPSignal40 = [derivativeBPSignal40 0];
% 
% % qrsLen = find(timeAxis == 132/sampFreq);
% 
% zeroPadDerBPSignal40  = [zeros(1,floor(windowLengthQRS/2)) derivativeBPSignal40 zeros(1,ceil(windowLengthQRS/2))];
% 
% nonLinearTransform = zeros(1,length(derivativeBPSignal40));
% for var1 = 1:1:length(derivativeBPSignal40)
%     nonLinearTransform(var1) = sum(abs(zeroPadDerBPSignal40(var1:var1 + windowLengthQRS))); % l(n)
% end
% 
% nonLinearTransform = nonLinearTransform - mean(nonLinearTransform);
% nonLinearTransform = nonLinearTransform/max(abs(nonLinearTransform));
% 
% adaptThresh = zeros(1,length(peaksLocationBandPass));
% for var1=1:1:length(peaksLocationBandPass)    
%     if var1 == 1
%         if peaksLocationBandPass(var1) - windowLengthQRS <= 0
%             adaptThresh(var1) = min(nonLinearTransform(1:peaksLocationBandPass(var1) + windowLengthQRS)) ... 
%                 + (max(nonLinearTransform(1:peaksLocationBandPass(var1) + windowLengthQRS)) - ... 
%                 min(nonLinearTransform(1:peaksLocationBandPass(var1) + windowLengthQRS)))/2.4;            
%         else
%             adaptThresh(var1) = min(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)) ... 
%                 + (max(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)) - ... 
%                 min(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)))/2.4;            
%         end        
%     elseif var1 == length(peaksLocationBandPass)
%         if peaksLocationBandPass(var1) + windowLengthQRS > length(nonLinearTransform)
%             adaptThresh(var1) = (adaptThresh(var1 - 1) + ... 
%                 min(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:end)) + ... 
%                 (max(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:end)) - ... 
%                 min(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:end)))/2.8)/2;            
%         else
%             adaptThresh(var1) = (adaptThresh(var1 - 1) + ... 
%                 min(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)) + ... 
%                 (max(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)) - ... 
%                 min(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)))/2.8)/2;
%         end
%     else
%         adaptThresh(var1) = (adaptThresh(var1 - 1) + ... 
%             min(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)) + ... 
%             (max(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)) - ... 
%             min(nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS)))/2.8)/2;        
%     end
% end
% 
% QRSTempOn = zeros(1,length(peaksLocationBandPass));
% QRSTempOff = zeros(1,length(peaksLocationBandPass));
% QRSTempOnLocation = zeros(1,length(peaksLocationBandPass));
% QRSTempOffLocation = zeros(1,length(peaksLocationBandPass));
% for var1=1:1:length(peaksLocationBandPass)
%     if var1 == 1
%         if peaksLocationBandPass(var1) - windowLengthQRS <= 0
%             thresholdIndex = nonLinearTransform(1:peaksLocationBandPass(var1) + windowLengthQRS) > adaptThresh(var1);
%             onFlag = 0;
%             offFlag = 1;
%             for var2=1:1:length(thresholdIndex)
%                 if thresholdIndex(var2) == 1 && offFlag == 1
%                     QRSTempOnLocation(var1) = var2;
%                     QRSTempOn(var1) = inSignal(QRSTempOnLocation(var1));
%                     onFlag = 1;
%                     offFlag = 0;
%                 elseif (thresholdIndex(var2) == 0 && onFlag == 1) || (onFlag == 1 && var2 == length(thresholdIndex))                    
%                     QRSTempOffLocation(var1) = var2;
%                     QRSTempOff(var1) = inSignal(QRSTempOffLocation(var1));
%                     onFlag = 0;
% %                     offFlag = 1;
%                 end
%             end
%         else
%             thresholdIndex = nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS) > adaptThresh(var1);
%             onFlag = 0;
%             offFlag = 1;
%             for var2=1:1:length(thresholdIndex)
%                 if thresholdIndex(var2) == 1 && offFlag == 1
%                     QRSTempOnLocation(var1) = peaksLocationBandPass(var1) - windowLengthQRS + var2 - 1;                    
%                     QRSTempOn(var1) = inSignal(QRSTempOnLocation(var1));
%                     onFlag = 1;
%                     offFlag = 0;
%                 elseif (thresholdIndex(var2) == 0 && onFlag == 1) || (onFlag == 1 && var2 == length(thresholdIndex))  
%                     QRSTempOffLocation(var1) = peaksLocationBandPass(var1) - windowLengthQRS + var2 - 1;
%                     QRSTempOff(var1) = inSignal(QRSTempOffLocation(var1));
%                     onFlag = 0;
% %                     offFlag = 1;
%                 end
%             end
%         end
%     elseif var1 == length(peaksLocationBandPass)        
%         if peaksLocationBandPass(var1) + windowLengthQRS > length(nonLinearTransform)
%             thresholdIndex = nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:end) > adaptThresh(var1);            
%             onFlag = 0;
%             offFlag = 1;
%             for var2=1:1:length(thresholdIndex)
%                 if thresholdIndex(var2) == 1 && offFlag == 1
%                     QRSTempOnLocation(var1) = peaksLocationBandPass(var1) - windowLengthQRS + var2 - 1;
%                     QRSTempOn(var1) = inSignal(QRSTempOnLocation(var1));
%                     onFlag = 1;
%                     offFlag = 0;
%                 elseif (thresholdIndex(var2) == 0 && onFlag == 1) || (onFlag == 1 && var2 == length(thresholdIndex))          
%                     QRSTempOffLocation(var1) = peaksLocationBandPass(var1) - windowLengthQRS + var2 - 1;
%                     QRSTempOff(var1) = inSignal(QRSTempOffLocation(var1));
%                     onFlag = 0;
% %                     offFlag = 1;
%                 end
%             end
%         else
%             thresholdIndex = nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS) > adaptThresh(var1);
%             onFlag = 0;
%             offFlag = 1;
%             for var2=1:1:length(thresholdIndex)
%                 if thresholdIndex(var2) == 1 && offFlag == 1
%                     QRSTempOnLocation(var1) = peaksLocationBandPass(var1) - windowLengthQRS + var2 - 1;
%                     QRSTempOn(var1) = inSignal(QRSTempOnLocation(var1));
%                     onFlag = 1;
%                     offFlag = 0;
%                 elseif (thresholdIndex(var2) == 0 && onFlag == 1) || (onFlag == 1 && var2 == length(thresholdIndex))        
%                     QRSTempOffLocation(var1) = peaksLocationBandPass(var1) - windowLengthQRS + var2 - 1;
%                     QRSTempOff(var1) = inSignal(QRSTempOffLocation(var1));
%                     onFlag = 0;
% %                     offFlag = 1;
%                 end
%             end            
%         end
%     else
%         thresholdIndex = nonLinearTransform(peaksLocationBandPass(var1) - windowLengthQRS:peaksLocationBandPass(var1) + windowLengthQRS) > adaptThresh(var1);
%         onFlag = 0;
%         offFlag = 1;
%         for var2=1:1:length(thresholdIndex)
%             if thresholdIndex(var2) == 1 && offFlag == 1
%                 QRSTempOnLocation(var1) = peaksLocationBandPass(var1) - windowLengthQRS + var2 - 1;
%                 QRSTempOn(var1) = inSignal(QRSTempOnLocation(var1));
%                 onFlag = 1;
%                 offFlag = 0;
%             elseif (thresholdIndex(var2) == 0 && onFlag == 1) || (onFlag == 1 && var2 == length(thresholdIndex))
%                 QRSTempOffLocation(var1) = peaksLocationBandPass(var1) - windowLengthQRS + var2 - 1;
%                 QRSTempOff(var1) = inSignal(QRSTempOffLocation(var1));
%                 onFlag = 0;
% %                 offFlag = 1;
%             end
%         end 
%     end
% end
% 
% 
% % figure
% % plot(inSignal,'b');
% % hold on;
% % plot(peaksLocationBandPass,peaksSignalBandPass,'r*');
% % hold on;
% % plot(QRSTempOnLocation,QRSTempOn,'g*');
% % hold on;
% % plot(QRSTempOffLocation,QRSTempOff,'g*');
% % hold on;
% 

% %% Locating Q, R and S
% 
% QPeakTemp = zeros(1,length(QRSTempOn));
% QPeakTempLocation = zeros(1,length(QRSTempOnLocation));
% RPeakTemp = zeros(1,length(peaksSignalBandPass));
% RPeakTempLocation = zeros(1,length(peaksLocationBandPass));
% SPeakTemp = zeros(1,length(QRSTempOff));
% SPeakTempLocation = zeros(1,length(QRSTempOffLocation));
% 
% for var1=1:1:length(peaksLocationBandPass)
%     [RPeakTemp(var1),RPeakTempLocation(var1)] = max(inSignal(QRSTempOnLocation(var1):QRSTempOffLocation(var1)));
%     RPeakTempLocation(var1) = QRSTempOnLocation(var1) + RPeakTempLocation(var1) - 1;
% 
%     [SPeakTemp(var1),SPeakTempLocation(var1)] = min(inSignal(RPeakTempLocation(var1):QRSTempOffLocation(var1)));
%     SPeakTempLocation(var1) = SPeakTempLocation(var1) + RPeakTempLocation(var1) - 1;
% 
%     [QPeakTemp(var1),QPeakTempLocation(var1)] = min(inSignal(QRSTempOnLocation(var1):RPeakTempLocation(var1)));
%     QPeakTempLocation(var1) = QPeakTempLocation(var1) + QRSTempOnLocation(var1) - 1;
% end
% 
% 
% QRSFeatures.QPeakTemp = QPeakTemp;
% QRSFeatures.QPeakTempLocation = QPeakTempLocation;
% QRSFeatures.RPeakTemp = RPeakTemp;
% QRSFeatures.RPeakTempLocation = RPeakTempLocation;
% QRSFeatures.SPeakTemp = SPeakTemp;
% QRSFeatures.SPeakTempLocation = SPeakTempLocation;
% 
% 
% %% Smoothing Filter filter
% 
% smoothingFilterSignal = filter([1/9 2/9 3/9 2/9 1/9],1,paddedInSignal);
% smoothingFilterSignal = smoothingFilterSignal(paddedWindowLength + 1:end);
% 
% % smoothingFilterSignal1 = smooth(paddedInSignal,windowLengthQRS + 1,'sgolay',5);
% % smoothingFilterSignal1 = smoothingFilterSignal1(paddedWindowLength + 1:end);
% % plot(inSignal,'b');
% % hold on;
% % plot(smoothingFilterSignal,'g');
% % hold on;
% % plot(smoothingFilterSignal1,'r');
% % hold on;
% 
% smoothingFilterSignal = smoothingFilterSignal - mean(smoothingFilterSignal);
% smoothingFilterSignal = smoothingFilterSignal/max(abs(smoothingFilterSignal));
% 
% derivativeSFSignal = diff(smoothingFilterSignal);
% derivativeSFSignal = [derivativeSFSignal 0];
% 
% derivativeSFSignal = derivativeSFSignal - mean(derivativeSFSignal);
% derivativeSFSignal = derivativeSFSignal/max(abs(derivativeSFSignal));
% 
% thresholdOn = zeros(1,length(QRSTempOnLocation));
% thresholdOff = zeros(1,length(QRSTempOffLocation));
% 
% thresholdOnLocation = zeros(1,length(QRSTempOnLocation));
% thresholdOffLocation = zeros(1,length(QRSTempOffLocation));
% QRSOnLocation = zeros(1,length(QRSTempOnLocation));
% QRSOffLocation = zeros(1,length(QRSTempOffLocation));
% for var1=1:1:length(QRSTempOnLocation)
%     firstPeak = QPeakTempLocation(var1);                               % check peaksLocationBandPass? Is this the first peak?
%     lastPeak = SPeakTempLocation(var1);                                % check peaksLocationBandPass? Is this the last peak?
%     [thresholdOn(var1),thresholdOnLocation(var1)] = max(abs(derivativeSFSignal(QRSTempOnLocation(var1):firstPeak)));
%     [thresholdOff(var1),thresholdOffLocation(var1)] = max(abs(derivativeSFSignal(lastPeak:QRSTempOffLocation(var1))));
%     
% %     thresholdOn(var1) = derivativeSFSignal(thresholdOnLocation(var1))/10;
% %     thresholdOff(var1) = derivativeSFSignal(thresholdOffLocation(var1))/10;
%     
%     thresholdOn(var1) = thresholdOn(var1)/10;
%     thresholdOff(var1) = thresholdOff(var1)/10;
%     
% %     figure;
% %     plot((QRSTempOnLocation(var1):QRSTempOffLocation(var1)),smoothingFilterSignal(QRSTempOnLocation(var1):QRSTempOffLocation(var1)),'k');
% %     hold on;
% %     plot((QRSTempOnLocation(var1):QRSTempOffLocation(var1)),derivativeSFSignal(QRSTempOnLocation(var1):QRSTempOffLocation(var1)),'b');
% %     hold on;
% % %     plot((QRSTempOnLocation(var1):QRSTempOffLocation(var1)),abs(derivativeSFSignal(QRSTempOnLocation(var1):QRSTempOffLocation(var1))),'g');
% % %     hold on;
% %     plot((QRSTempOnLocation(var1):firstPeak),derivativeSFSignal(QRSTempOnLocation(var1):firstPeak),'r');
% %     hold on;
% %     plot((lastPeak:QRSTempOffLocation(var1)),derivativeSFSignal(lastPeak:QRSTempOffLocation(var1)),'r');
% %     hold on;
% %     plot((QRSTempOnLocation(var1):QRSTempOffLocation(var1)),thresholdOn(var1)*ones(1,length(derivativeSFSignal(QRSTempOnLocation(var1):QRSTempOffLocation(var1)))),'g');
% %     hold on;
% %     plot((QRSTempOnLocation(var1):QRSTempOffLocation(var1)),thresholdOff(var1)*ones(1,length(derivativeSFSignal(QRSTempOnLocation(var1):QRSTempOffLocation(var1)))),'b');
% %     hold on;
%     
%     
%     thresholdOnLocation(var1) = thresholdOnLocation(var1) + QRSTempOnLocation(var1) - 1;
%     thresholdOffLocation(var1) = thresholdOffLocation(var1) + lastPeak - 1;
%     
%     searchFlag1 = 1;
%     currentLocation = firstPeak;
%     while searchFlag1
%         if derivativeSFSignal(currentLocation) < thresholdOn(var1) || currentLocation == QRSTempOnLocation(var1)
%             QRSOnLocation(var1) = currentLocation;
%             searchFlag1 = 0;
%         else
%             currentLocation = currentLocation - 1;
%         end
%             
%     end
%     
%     searchFlag2 = 1;
%     currentLocation = QRSTempOffLocation(var1);
%     while searchFlag2
%         if derivativeSFSignal(currentLocation) > thresholdOff(var1) || currentLocation == lastPeak
%             if currentLocation == lastPeak
%                 QRSOffLocation(var1) = QRSTempOffLocation(var1);
%                 searchFlag2 = 0;
%             else
%                 QRSOffLocation(var1) = currentLocation;
%                 searchFlag2 = 0;
%             end
%         else
%             currentLocation = currentLocation - 1;
%         end
%             
%     end
%     
%     
%     
% %     searchFlag1 = 1;
% % %     currentLocation = thresholdOnLocation(var1);
% %     currentLocation = QPeakTempLocation(var1);
% %     while searchFlag1
% % %         if derivativeSFSignal(currentLocation) < thresholdOn(var1) || currentLocation == 1
% %         if derivativeSFSignal(currentLocation) < thresholdOn(var1) || currentLocation == QRSTempOnLocation(var1)
% %             QRSOnLocation(var1) = currentLocation;
% %             searchFlag1 = 0;
% %         else
% %             currentLocation = currentLocation - 1;
% %         end
% %             
% %     end
% %     
% %     searchFlag2 = 1;
% % %     currentLocation = thresholdOffLocation(var1);
% %     currentLocation = QRSTempOffLocation(var1);
% %     while searchFlag2
% % %         if derivativeSFSignal(currentLocation) > thresholdOff(var1) || currentLocation == length(inSignal)
% %         if derivativeSFSignal(currentLocation) > thresholdOff(var1) || currentLocation == QRSTempOffLocation
% %             QRSOffLocation(var1) = currentLocation;
% %             searchFlag2 = 0;
% %         else
% %             currentLocation = currentLocation + 1;
% %         end
% %             
% %     end
% end
% 
% 
% QRSOn = inSignal(QRSOnLocation);
% QRSOff = inSignal(QRSOffLocation);
% 
% % QRSOn = smoothingFilterSignal(QRSOnLocation);
% % QRSOff = smoothingFilterSignal(QRSOffLocation);
% 
% 
% QRSFeatures.QRSOn = QRSOn;
% QRSFeatures.QRSOff = QRSOff;
% 
% % figure
% % plot(inSignal,'b');
% % hold on;
% % plot(peaksLocationBandPass,peaksSignalBandPass,'r*');
% % hold on;
% % plot(QRSOnLocation,QRSOn,'g*');
% % hold on;
% % plot(QRSOffLocation,QRSOff,'g*');
% % hold on;
% 
