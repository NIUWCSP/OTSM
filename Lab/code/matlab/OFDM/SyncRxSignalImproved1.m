function startIdx = SyncRxSignalImproved1(rxFrame, overSampFactor, numFFT,M_mod,N,M)
%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix

%% Definitions 同步找開頭結尾
numShortPreambleSamples = 16     * overSampFactor;
numLongPreambleSamples  = numFFT * overSampFactor;

thresholdCoarse = 0.6;%改回原值
thresholdFine   = 0.6;%改回原值

frameLen = length(rxFrame);
%% Set start index to an invalid number
startIdx = -1;
% M_mod: size of QAM constellation
M_bits = log2(M_mod);

%% Construct the syncSig to be used for fine tuning

%%調變同步資料
SyncBits = GetSyncBits();%確保正確解讀接收到的數據
SyncSymb_tilda=QamAndTilda(SyncBits,M_mod,M_bits,N,M,Wn);

syncSig = SyncSymb_tilda;
%syncSig = ifft(SyncSymb_tilda) * sqrt(length(SyncSymb_tilda));

%% Cross correlate different segments of the Rx signal, and the sync signal
corrShortCoarse  = zeros(1, frameLen);
corrFine  = zeros(1, frameLen);
% Region of interest, for visualization and plotting
roi = zeros(1, frameLen);

for i = 1: frameLen - 2*numShortPreambleSamples - 3*numLongPreambleSamples
     % Grab A1
     initIdx = i;
     seg1 = rxFrame( initIdx : initIdx + numShortPreambleSamples - 1 );
     % Grab first 32 symbols of B
     initIdx = initIdx + 2*numShortPreambleSamples;
     seg2 = rxFrame( initIdx : initIdx + numShortPreambleSamples - 1 );
%      seg1=zeros(0);    
%      for j = 0:numShortPreambleSamples/size(QamSyncBits,1)-1
%         seg1 = [seg1
%                 rxFrame( initIdx+j*M : initIdx+j*M + size(QamSyncBits,2) - 1 ,1)];
%      end
%      % Grab first 32 symbols of B (8+8 中間空64格)
%      initIdx = initIdx + 2*numShortPreambleSamples;
%      seg2=zeros(0);
%      for j = 0:numShortPreambleSamples/size(QamSyncBits,1)-1
%         seg2 = [seg2
%                 rxFrame( initIdx+j*M : initIdx+j*M + size(QamSyncBits,2) - 1 ,1)];
%      end
     % Normalization factors
     seg1Avg = sqrt(sum(abs(seg1).^2));
     seg2Avg = sqrt(sum(abs(seg2).^2));
     % Cross correlation of short preambles
     corrShortCoarse(i) = abs(sum(seg1 .* conj(seg2))) / (seg1Avg * seg2Avg);
     corrShortCoarse1(i) = abs(sum(seg1 .* conj(seg2)));
     corrShortCoarse2(i) = (seg1Avg * seg2Avg);
     % If short preambles are correlated, we are
     % at a potential coarse begining of an OTSM frame, therefore check
     % Fine synchronization criterion
    if (corrShortCoarse(i) > thresholdCoarse)
         % Grab B pilot
%         initIdx  = i + numFFT*2+40+size(QamSyncBits,2);
%         segPilot=zeros(0);
%         for k = 0:numLongPreambleSamples/size(QamSyncBits,1)-1
%             segPilot = [segPilot
%                         rxFrame( initIdx+k*M+1 : initIdx+k*M + size(QamSyncBits,2) ,1)];
%         end
         
          initIdx  = i + 2*numShortPreambleSamples;
          segPilot = rxFrame(initIdx : initIdx + numLongPreambleSamples - 1);

         % Normalization factors
         segPilotAvg = sqrt(sum(abs(segPilot).^2));
         syncSigAvg  = sqrt(sum(abs(syncSig) .^2));
         corrFine(i) = abs(sum(segPilot .* conj(syncSig))) / ...
                       (segPilotAvg * syncSigAvg);
         if corrFine(i) > thresholdFine
             % Mark this index in the region of interest
             roi(i) = 1;            
             if startIdx == -1
                startIdx = i;
             end
         end
     end
end

if (startIdx == -1)
    disp('No OTSM frame was found')
else
    disp(['OTSM frame startIdx = ', num2str(startIdx)]);
end

subplot(231);
set(gcf, 'Position', [1300, 400, 600, 500]);
tAxis = (0:frameLen-1)/((1.4e6 * overSampFactor)/10^3);
hold on
plot(tAxis, abs(corrShortCoarse))
plot(tAxis, abs(corrFine))
plot(tAxis, roi)
title('sync plot');
xlabel('time [ms]')
legend('corrShortCoarse', 'corrFine', 'ROI')
hold off
