function RxDataBits = Receiver(RxSignal)

% Setting parameters
NumFFT = 128; %V3
NumSyncPreamble = 32; %V3
NumCP = 16; %V3

NumDataOfdmSymb = 18;
NumDataSubcarrier = 108;
RxSignalExt(:,1)=RxSignal;
figure(2);clf;

%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
% M: number of subcarriers in frequency
M = 64;
% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame

%% Receiver
NumSyncSymb = NumSyncPreamble*2 + NumFFT;
NumPilotSymb = NumFFT * 2;
NumDataSymb = (NumFFT+NumCP) * NumDataOfdmSymb;
NumRadioFrame = NumSyncSymb + NumPilotSymb + NumDataSymb;

StartIdx = SyncRxSignalImproved1(RxSignalExt, 1, NumFFT);

%%JF加入重新賦值=1避免StartIdx == -1時直接中斷程式：
global NoFoundDataTimes;
if(StartIdx == -1)
    StartIdx = 1;
    NoFoundDataTimes=NoFoundDataTimes+1;
elseif(StartIdx+NumRadioFrame-1 >= 31680)
    StartIdx = 1;
    NoFoundDataTimes=NoFoundDataTimes+1;
end
RxSignalRadioFrame = RxSignalExt(StartIdx:StartIdx+NumRadioFrame-1);

% Pilot OFDM symbol
PilotOfdmSymb = reshape(RxSignalRadioFrame(NumSyncSymb+1:NumSyncSymb+NumPilotSymb), [], 2);

 % 3GPP channel model
 max_speed=500;  % km/hr
[chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);

 %% channel output%%%%% 
[G,gs]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);

%% Demodulation
% Estimate carrier frequency offset
RxPilotSymb = OfdmSignalDemodulation(PilotOfdmSymb, NumFFT, 0, NumDataSubcarrier);
XCorrPilot = RxPilotSymb(:,1)' * RxPilotSymb(:,2);
EpsEst = 1/(2*pi) * atan(imag(XCorrPilot)/real(XCorrPilot));

% Estimate carrier freqnecy offset
RxSigalRadioFrameCmpCFO = RxSignalRadioFrame .* ...
    exp(-1j*2*pi*EpsEst/NumFFT * (0:length(RxSignalRadioFrame)-1)');

% Reobtain pilot data
PilotOfdmSymb = reshape( ...
    RxSigalRadioFrameCmpCFO(NumSyncSymb+1:NumSyncSymb+NumPilotSymb), [], 2);
% Data OFDM symbol
DataOfdmSymb = reshape( ...
    RxSigalRadioFrameCmpCFO(NumSyncSymb+NumPilotSymb+1:end), ...
    [], NumDataOfdmSymb);

RxPilotSymb = OfdmSignalDemodulation(PilotOfdmSymb, NumFFT, 0, NumDataSubcarrier);
RxDataSymb = OfdmSignalDemodulation(DataOfdmSymb, NumFFT, NumCP, NumDataSubcarrier);

% Channel estimation and equalization
BpskModObj = comm.BPSKModulator('PhaseOffset', pi/4);
PilotBits = GetPilotBits();
TxPilotSymb = step(BpskModObj,PilotBits);
ChanEst = RxPilotSymb(:,1) ./ TxPilotSymb;
global RxDataSymbEq;
RxDataSymbEq = RxDataSymb ./ repmat(ChanEst, 1, NumDataOfdmSymb);
subplot(232);plot(10*log10(abs(ChanEst).^2)-min(10*log10(abs(ChanEst).^2)));title('channel estimation');
subplot(233);plot(RxDataSymb(:),'*');axis equal;title('');title('scatter before equalization');axis square;
subplot(234);plot(RxDataSymbEq(:).*exp(-1i*pi/4),'.');axis([-1.5,1.5,-1.5,1.5]);title('scatter after equalization'); axis square;
% Demodulation
BpskDemodulatorObj = comm.BPSKDemodulator( ...
    'PhaseOffset', pi/4, ...
    'DecisionMethod', 'Hard decision');
global RxDataBits;
RxDataBits = zeros(size(RxDataSymbEq));
for idx = 1:NumDataOfdmSymb
    RxDataBits(:,idx) = step(BpskDemodulatorObj,RxDataSymbEq(:,idx));
end
