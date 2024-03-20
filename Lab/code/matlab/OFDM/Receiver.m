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
M_data = M-delay_spread;%64-24=40
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));%64*40=2,560
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;



% Time and frequency resources
car_fre = 2.4*10^9;% Carrier frequency
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
%%  Receiver
RxSignalExt(:,1)=RxSignal;
%PilotNumDataSubcarrier =64;
%DataNumDataSubcarrier =64;

        NumSyncSymb =  NumSyncPreamble*2+size(GetSyncBits,2)/2;%
        NumPilotSymb = 0;%PilotBits為128個Bits QAM後會除2
        NumDataSymb = N*M;
        NumRadioSymb = NumSyncSymb + NumPilotSymb + NumCP + NumDataSymb;
        
        figure(2);clf;
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
        RxSignalRadioFrame = RxSignalExt(StartIdx + NumSyncSymb + NumPilotSymb + NumCP:StartIdx+NumRadioSymb-1);

        %RxSignalDataFrame=OtsmSignalDemodulation(RxSignalRadioFrame, NumFFT, 0, DataNumDataSubcarrier,M_mod);
        %RxSignalDataFrame=reshape(RxSignalDataFrame,[],1);

         %% OTSM Reobtain Pilot%%%%
                
         %Pilot side
         RxSignalGrid=reshape(RxSignalRadioFrame,N,M);
         Y_OTSM_Pilot=reshape(RxSignalGrid(M_data+sqrt(size(GetPilotBits,2)/2)+1:M_data+8*2,1:sqrt(size(GetPilotBits,2)/2)),[],1);%傳送資料時Pilot和Sync皆以8*8的形式傳送

         %%把pilot的8*8資料 拆成4*8、4*8
         Y_OTSM_PilotSymb=zeros(size(Y_OTSM_Pilot,1)/2,2);
         for i=0:sqrt(size(Y_OTSM_Pilot,1))-1
               Y_OTSM_PilotSymb(i*4+1:i*4+4,1) = Y_OTSM_Pilot(i*8+1:i*8+4,1);
               Y_OTSM_PilotSymb(i*4+1:i*4+4,2) = Y_OTSM_Pilot(i*8+5:i*8+8,1);
         end

        % Estimate carrier frequency offset
        [RxDataSymbEq,ChanEst] = channel_est(N,M,M_mod,NumFFT,RxSignalRadioFrame,Y_OTSM_PilotSymb);

        %EQ測試用
                RxSymbEq=reshape([RxDataSymbEq;
                                  zeros(N-M_data,M)],[],1);


         %% OTFS channel generation%%%%
         % 3GPP channel model
         max_speed=500;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        
         %% channel output%%%%% 
        [G,gs]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);

        r=zeros(N*M,1);
        noise = sqrt(sigma_2(iesn0)/2) * (randn(size(RxSignalRadioFrame)) + 1i*randn(size(RxSignalRadioFrame)));
        l_max=max(delay_taps);
        for q=0:N*M-1
            for l=0:l_max
                if(q>=l)
                    r(q+1)=r(q+1)+gs(l+1,q+1)*RxSymbEq(q-l+1);  %equation (24) in [R1]
                end
            end
        end
        r=r+noise;

         

        %% OTSM demodulation%%%%

            %主要解調變
                Y_tilda=reshape(r,M,N);     %equation (11) in [R1]
                Y = Y_tilda*Wn;             %equation (12) in [R1]
            
         
        
        %% Generate the block-wise channel matrices in the delay-time and the time-frequency domain
        [Gn_block_matrix,Tn_block_matrix,zn_block_vector,H_t_f]=Generate_Matched_Filter_GS_matrices(N,M,G,r,ChanEst);
        
         %% GS SOR Iterative detection
        
        n_ite_MRC=50; % maximum number of detector iterations
        omega=1; %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        if(M_mod==64)
            omega=0.25;
        end
        decision=1; %1-hard decision, 0-soft decision
        [est_info_bits_MFGS,det_iters_MFGS,~] = Matched_Filter_GS_detector(N,M,M_mod,sigma_2(iesn0),data_grid,Y,H_t_f,n_ite_MRC,omega,Tn_block_matrix,Gn_block_matrix,zn_block_vector,r,Wn,decision);
        [est_info_bits_1tap,~] = TF_single_tap_equalizer(N,M,M_mod,sigma_2(iesn0),data_grid,Y,H_t_f,Wn);
        [est_info_bits_LMMSE,~] = Block_LMMSE_detector(N,M,M_mod,sigma_2(iesn0),data_grid,Gn_block_matrix,r,Wn);
        RxDataBits=0;
       
