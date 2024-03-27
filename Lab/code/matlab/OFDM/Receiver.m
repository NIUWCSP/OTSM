function [RxDataBits,est_info_bits_MFGS,det_iters_MFGS,est_info_bits_1tap,est_info_bits_LMMSE] = Receiver(RxSignal, sigma, N, M, M_mod)


%%廣域變數宣告
% global iesn0

%% OTFS parameters%%%%%%%%%%

% average energy per data symbol
% eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));


%% delay-Doppler grid symbol placement
% max delay spread in the channel
delay_spread = M/(8/3);%40*64是資料部分 剩下是Pilot跟Sync
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-delay_spread;%64-24=40
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;


% Time and frequency resources
% car_fre = 2.4*10^9;% Carrier frequency 原先4*10^9
% delta_f = 15*10^3; % subcarrier spacing: 15 KHz
% T = 1/delta_f; %one time symbol duration in OTFS frame




%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix

% Setting parameters

NumSyncPreamble = 32; %V3 同步的前綴，Preamble：防干擾+同步+通道估測(已知的頻域資料)
NumCP = 16; %V3 CP：循環前綴(NumFFT = 128後16貼回前面)，CP：避免ISI(多路徑干擾)(未知的時域訊號)
%%  Receiver
RxSignalExt(:,1)=RxSignal;

        NumSyncSymb =  NumSyncPreamble*2+128;%128是調變後的Sync值
        NumDataSymb = N*M;
        NumRadioSymb = NumSyncSymb + NumCP + NumDataSymb;
        
        %%畫圖準備
        figure(3);clf;
        
        %StartIdx = SyncRxSignalImproved2(RxSignalExt,M_mod,N,M);
        StartIdx = SyncRxSignalImproved1(RxSignalExt, 1 ,M_mod,N,M,Wn);

        %%JF加入重新賦值=1避免StartIdx == -1時直接中斷程式：
        global NoFoundDataTimes;
        if(StartIdx == -1)
            StartIdx = 1;
            NoFoundDataTimes=NoFoundDataTimes+1;
        elseif(StartIdx+NumDataSymb-1 >= length(RxSignalExt))
            StartIdx = 1;
            NoFoundDataTimes=NoFoundDataTimes+1;
        end
        RxSignalRadioFrame = RxSignalExt(StartIdx + NumSyncSymb  + NumCP:StartIdx+NumRadioSymb-1);

         %% OTSM Reobtain Pilot%%%%
                
         %Pilot side
         RxSignalGrid=reshape(RxSignalRadioFrame,N,M);
         Y_OTSM_Pilot=zeros(N,2);
         for i=1:sqrt(size(Y_OTSM_Pilot,1))
         %Y_OTSM_Pilot=reshape(RxSignalGrid(M_data+sqrt(size(GetPilotBits,2)/2)+1:M_data+8*2,1:sqrt(size(GetPilotBits,2)/2)),[],1);%傳送資料時Pilot和Sync皆以8*8的形式傳送
         Y_OTSM_Pilot((i-1)*8+1:(i-1)*8+8,1)=RxSignalGrid(M_data+1:M_data+sqrt(size(GetPilotBits,2)/2),i);
         Y_OTSM_Pilot((i-1)*8+1:(i-1)*8+8,2)=RxSignalGrid(M_data+sqrt(size(GetPilotBits,2)/2)*2+1:M_data+sqrt(size(GetPilotBits,2)/2)*3,i);
         end
         %%把pilot的8*8資料 拆成4*8、4*8

%             %% OTFS channel generation%%%%  for math model
%          % 3GPP channel model
%          max_speed=500;  % km/hr
%         [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
%         
%          %% channel output%%%%% 
%         [G,gs,l_max]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);
% 
%        
% 
% 
% 
% 
%         r=zeros(N*M,1);
%         noise = sigma * (randn(size(RxSignalRadioFrame)) + 1i*randn(size(RxSignalRadioFrame)));
%         l_max=max(delay_taps);
%         for q=0:N*M-1
%             for l=0:l_max
%                 if(q>=l)
%                     r(q+1)=r(q+1)+gs(l+1,q+1)*RxSymbEq(q-l+1);  %equation (24) in [R1]
%                 end
%             end
%         end
%         r=r+noise;
% 
%          %%% Estimating the practical channel model 
%          %%% outputs [G,gs,l_max]
 % Estimate carrier frequency offset
        [RxDataSymbEq,RxSigalRadioFrameCmpCFO,G] = channel_est(N,M,M_mod,RxSignalRadioFrame,Y_OTSM_Pilot,0);

        %EQ測試用
                %RxSymbEq=reshape([RxDataSymbEq;
                                  %zeros(N-M_data,M)],[],1);
        r = reshape(RxSigalRadioFrameCmpCFO,[],1);
        %% OTSM demodulation%%%%

            %主要解調變
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
        [est_info_bits_MFGS,det_iters_MFGS,~] = Matched_Filter_GS_detector(N,M,M_mod,sigma,data_grid,Y,H_t_f,n_ite_MRC,omega,Tn_block_matrix,Gn_block_matrix,zn_block_vector,r,Wn,decision);
        [est_info_bits_1tap,~] = TF_single_tap_equalizer(N,M,M_mod,sigma,data_grid,Y,H_t_f,Wn);
        [est_info_bits_LMMSE,~] = Block_LMMSE_detector(N,M,M_mod,sigma,data_grid,Gn_block_matrix,r,Wn);
        RxDataBits=0;
       
