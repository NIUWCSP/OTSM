function [RxDataBits,est_info_bits_MFGS,det_iters_MFGS,est_info_bits_1tap,est_info_bits_LMMSE] = Receiver(RxSignal)

%%廣域變數宣告
global iesn0

%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));


%% delay-Doppler grid symbol placement
% max delay spread in the channel
delay_spread = M/(8/3);%40*64是資料部分 剩下是Pilot跟Sync
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-delay_spread;
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;



% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 10:2.5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;



%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix

% Setting parameters
NumFFT = 64; %V3  FFT轉換的點數
NumSyncPreamble = 32; %V3 同步的前綴，Preamble：防干擾+同步+通道估測(已知的頻域資料)
NumCP = 16; %V3 CP：循環前綴(NumFFT = 128後16貼回前面)，CP：避免ISI(多路徑干擾)(未知的時域訊號)

PilotSymb=size(GetPilotBits,2)/2;%PilotBits為128個Bits QAM後會除2
RxSignalExt(:,1)=RxSignal;
%PilotNumDataSubcarrier =64;
%DataNumDataSubcarrier =64;

        %% Receiver 

        NumDataSymb = N*M;
        
        %StartIdx = SyncRxSignalImproved2(RxSignalExt,M_mod,N,M);
        StartIdx = SyncRxSignalImproved1(RxSignalExt, 1, NumFFT,M_mod,N,M);

        %%JF加入重新賦值=1避免StartIdx == -1時直接中斷程式：
        global NoFoundDataTimes;
        if(StartIdx == -1)
            StartIdx = 1;
            NoFoundDataTimes=NoFoundDataTimes+1;
        elseif(StartIdx+NumDataSymb-1 >= length(RxSignalExt))
            StartIdx = 1;
            NoFoundDataTimes=NoFoundDataTimes+1;
        end
        RxSignalRadioFrame = RxSignalExt(StartIdx:StartIdx+NumDataSymb-1);

        %RxSignalDataFrame=OtsmSignalDemodulation(RxSignalRadioFrame, NumFFT, 0, DataNumDataSubcarrier,M_mod);
        %RxSignalDataFrame=reshape(RxSignalDataFrame,[],1);

            

         %% OTFS channel generation%%%%
         % 3GPP channel model
         max_speed=500;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        
         %% channel output%%%%% 
        [G,~]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);

         %% OTSM Reobtain Pilot%%%%
                
         %Pilot side
         RxSignalGrid=reshape(RxSignalRadioFrame,N,M);
         Y_OTSM_Pilot=RxSignalGrid(M_data+sqrt(PilotSymb)+1:M_data+sqrt(PilotSymb)*2,1:sqrt(PilotSymb));%傳送資料時Pilot和Sync皆以8*8的形式傳送

        % Estimate carrier frequency offset
XCorrPilot = Y_OTSM_Pilot(:,1)' * Y_OTSM_Pilot(:,2);%導頻符號的交叉相關
EpsEst = 1/(2*pi) * atan(imag(XCorrPilot)/real(XCorrPilot));%頻率偏移的估計

% Estimate carrier freqnecy offset
RxSigalRadioFrameCmpCFO = RxSignalRadioFrame .* ...
    exp(-1j*2*pi*EpsEst/NumFFT * (0:length(RxSignalRadioFrame)-1)');%接收訊號進行CFO校正
RxSignalRadioGridCFO = reshape(RxSigalRadioFrameCmpCFO,N,M);

% Reobtain pilot data
PilotOtsmSymb = RxSignalRadioGridCFO(M_data+sqrt(PilotSymb)+1:M_data+sqrt(PilotSymb)*2,1:sqrt(PilotSymb));%將所選的一段導頻資料重新組織成一個矩陣，其中每列有兩個元素

% Data OFDM symbol
DataOtsmSymb = RxSignalRadioGridCFO(1:M_data,1:M);


% Channel estimation and equalization
PilotBits = GetPilotBits();
Tx_PilotSymb = qammod(reshape(PilotBits,M_bits,[]),M_mod,'gray','InputType','bit');%產生Tx的 BPSK符號序列
ChanEst = reshape(PilotOtsmSymb,1,[]) ./ Tx_PilotSymb;%通道估计
global RxDataSymbEq;
RxDataSymbEq = DataOtsmSymb ./ repmat(ChanEst, M_data,1);
subplot(232);plot(10*log10(abs(ChanEst).^2)-min(10*log10(abs(ChanEst).^2)));title('channel estimation');%繪製通道估計的幅度譜
subplot(233);plot(DataOtsmSymb(:),'*');axis equal;title('scatter before equalization');axis square;
subplot(234);plot(RxDataSymbEq(:).*exp(-1i*pi/4),'.');axis([-1.5,1.5,-1.5,1.5]);title('scatter after equalization'); axis square;%*exp(-1i*pi/4) 的作用是進行相位調整

        %% OTSM demodulation%%%%
                RxSymbEq=[RxDataSymbEq;
                          zeros(N-M_data,M)  ];
                Y_tilda=reshape(RxSymbEq,M,N);     %equation (11) in [R1]
                Y = Y_tilda*Wn;             %equation (12) in [R1]
               
         
        
        %% Generate the block-wise channel matrices in the delay-time and the time-frequency domain
        [Gn_block_matrix,Tn_block_matrix,zn_block_vector,H_t_f]=Generate_Matched_Filter_GS_matrices(N,M,G,RxSignalRadioFrame);
        
         %% GS SOR Iterative detection
        
        n_ite_MRC=50; % maximum number of detector iterations
        omega=1; %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        if(M_mod==64)
            omega=0.25;
        end
        decision=1; %1-hard decision, 0-soft decision
        [est_info_bits_MFGS,det_iters_MFGS,data_MFGS] = Matched_Filter_GS_detector(N,M,M_mod,sigma_2(iesn0),data_grid,Y,H_t_f,n_ite_MRC,omega,Tn_block_matrix,Gn_block_matrix,zn_block_vector,RxSignalRadioFrame,Wn,decision);
        [est_info_bits_1tap,data_1tap] = TF_single_tap_equalizer(N,M,M_mod,sigma_2(iesn0),data_grid,Y,H_t_f,Wn);
        [est_info_bits_LMMSE,data_LMMSE] = Block_LMMSE_detector(N,M,M_mod,sigma_2(iesn0),data_grid,Gn_block_matrix,RxSignalRadioFrame,Wn);
        RxDataBits=0;
       
