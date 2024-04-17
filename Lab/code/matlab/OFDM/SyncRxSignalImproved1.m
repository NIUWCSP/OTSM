function startIdx = SyncRxSignalImproved1(rxFrame, overSampFactor,M_mod,N,M,Wn)

%% Definitions 同步找開頭結尾
numShortPreambleSamples = 32     * overSampFactor;
numLongPreambleSamples  = M*2 * overSampFactor;

thresholdCoarse = 0.9;%改回原值
thresholdFine   = 0.6;%改回原值

frameLen = length(rxFrame);
%% Set start index to an invalid number
startIdx = -1;
% M_mod: size of QAM constellation
M_bits = log2(M_mod);

%% Construct the syncSig to be used for fine tuning

%%調變同步資料

SyncBits = GetSyncBits();%確保正確解讀接收到的數據
QamSyncBits=reshape(qammod(reshape(SyncBits,M_bits,size(SyncBits,2)/M_bits), M_mod,'gray','InputType','bit'),[],1);
QamSync_ISI = CompeteISI(QamSyncBits,0); %128*1

syncSig = QamSync_ISI;
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
ylabel('autocorrelation')
legend('corrShortCoarse', 'corrFine', 'ROI')
hold off
