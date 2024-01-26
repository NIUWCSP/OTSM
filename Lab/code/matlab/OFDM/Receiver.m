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
% M_mod: size of QAM constellation
M_mod = 4;

%% delay-Doppler grid symbol placement
% max delay spread in the channel
delay_spread = M/16;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-delay_spread;
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;


% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame

%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix

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

 %% OTFS channel generation%%%%
 % 3GPP channel model
 max_speed=500;  % km/hr
[chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);

 %% channel output%%%%% 
[G,gs]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);

 r=zeros(N*M,1);
 l_max=max(delay_taps);
         for q=0:N*M-1
            for l=0:l_max
                if(q>=l)
                    r(q+1)=r(q+1)+gs(l+1,q+1)*RxSignal(q-l+1);  %equation (24) in [R1]
                end
            end
        end
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

%% OTSM demodulation%%%%
        Y_tilda=reshape(r,M,N);     %equation (11) in [R1]
        Y = Y_tilda*Wn;             %equation (12) in [R1]

%% Generate the block-wise channel matrices in the delay-time and the time-frequency domain
[Gn_block_matrix,Tn_block_matrix,zn_block_vector,H_t_f]=Generate_Matched_Filter_GS_matrices(N,M,G,r);

 %% GS SOR Iterative detection
        
        n_ite_MRC=50; % maximum number of detector iterations
        omega=1; %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        if(M_mod==64)
            omega=0.25;
        end
        decision=1; %1-hard decision, 0-soft decision
        [est_info_bits_MFGS,det_iters_MFGS,data_MFGS] = Matched_Filter_GS_detector(N,M,M_mod,0,data_grid,Y,H_t_f,n_ite_MRC,omega,Tn_block_matrix,Gn_block_matrix,zn_block_vector,r,Wn,decision);
        [est_info_bits_1tap,data_1tap] = TF_single_tap_equalizer(N,M,M_mod,0,data_grid,Y,H_t_f,Wn);
        [est_info_bits_LMMSE,data_LMMSE] = Block_LMMSE_detector(N,M,M_mod,0,data_grid,Gn_block_matrix,r,Wn);

         %% errors count%%%%%
        errors_MFGS = sum(xor(est_info_bits_MFGS,trans_info_bit));
        errors_1tap = sum(xor(est_info_bits_1tap,trans_info_bit));
        errors_LMMSE = sum(xor(est_info_bits_LMMSE,trans_info_bit));
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
