function [est_info_bits_MFGS,det_iters_MFGS,est_info_bits_1tap,est_info_bits_LMMSE] = Receiver(RxSignal, sigma, N, M, M_mod,M_bits)


%%廣域變數宣告
% global iesn0

%% OTFS parameters%%%%%%%%%%

%% delay-Doppler grid symbol placement
% max delay spread in the channel
delay_spread = M/(8/3);%40*64是資料部分 剩下是Pilot跟Sync
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-delay_spread;%64-24=40
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
DelayPilotSymb=sqrt(size(GetPilotBits,2)/2);


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
% %%把兩個Pilot分別放入2行
%          Y_OTSM_Pilot=zeros(N,2);
%          for i=1:sqrt(size(Y_OTSM_Pilot,1))
%          %Y_OTSM_Pilot=reshape(RxSignalGrid(M_data+DelayPilotSymb+1:M_data+8*2,1:DelayPilotSymb),[],1);%傳送資料時Pilot和Sync皆以8*8的形式傳送
%          Y_OTSM_Pilot((i-1)*8+1:(i-1)*8+8,1)=RxSignalGrid(M_data+1:M_data+DelayPilotSymb,i);
%          Y_OTSM_Pilot((i-1)*8+1:(i-1)*8+8,2)=RxSignalGrid(M_data+DelayPilotSymb*2+1:M_data+DelayPilotSymb*3,i);
%          end

%%把兩個Pilot切割成4*n後交叉放入2行
         Y_OTSM_Pilot=zeros(N,2);
         for i=1:sqrt(size(Y_OTSM_Pilot,1))
         %Y_OTSM_Pilot=reshape(RxSignalGrid(M_data+DelayPilotSymb+1:M_data+8*2,1:DelayPilotSymb),[],1);%傳送資料時Pilot和Sync皆以8*8的形式傳送
         Y_OTSM_Pilot((i-1)*4+1:(i-1)*4+4,1)=RxSignalGrid(M_data+1:M_data+4,i);
         Y_OTSM_Pilot((i-1)*4+1:(i-1)*4+4,2)=RxSignalGrid(M_data+5:M_data+8,i);
         end
         for i=1:sqrt(size(Y_OTSM_Pilot,1))
         %Y_OTSM_Pilot=reshape(RxSignalGrid(M_data+DelayPilotSymb+1:M_data+8*2,1:DelayPilotSymb),[],1);%傳送資料時Pilot和Sync皆以8*8的形式傳送
         Y_OTSM_Pilot((i+7)*4+1:(i+7)*4+4,1)=RxSignalGrid(M_data+DelayPilotSymb+1:M_data+DelayPilotSymb+4,i);
         Y_OTSM_Pilot((i+7)*4+1:(i+7)*4+4,2)=RxSignalGrid(M_data+DelayPilotSymb+5:M_data+DelayPilotSymb+8,i);
         end
        
         %把pilot的8*8資料 拆成4*8、4*8

 % Estimate carrier frequency offset
        [RxDataSymbEq,RxSignalRadioFrameCmpCFO,G] = channel_est(N,M,M_mod,RxSignalRadioFrame,Y_OTSM_Pilot,0);
        r = reshape(RxSignalRadioFrameCmpCFO,[],1);
        
        %% OTSM demodulation%%%%

            %主要解調變
                Y_tilda=reshape(r,M,N);     %equation (11) in [R1]
                Y = Y_tilda*Wn;             %equation (12) in [R1]
            
         
        
        %% Generate the block-wise channel matrices in the delay-time and the time-frequency domain
        [Gn_block_matrix,Tn_block_matrix,zn_block_vector,H_t_f]=Generate_Matched_Filter_GS_matrices(N,M,G,r);
        %%偵錯
        if (isnan(sum(zn_block_vector(:))))
            est_info_bits_MFGS=zeros(M_data*M*M_bits,1);
            est_info_bits_1tap=zeros(M_data*M*M_bits,1);
            est_info_bits_LMMSE=zeros(M_data*M*M_bits,1);
            det_iters_MFGS=0;
            return;
        end
         %% GS SOR Iterative detection
        
        n_ite_MRC=50; % maximum number of detector iterations
        omega=1; %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        if(M_mod==64)
            omega=0.25;
        end
        figure(5);
        set(gcf,'name','三種偵測器EQ圖','Position', [650 50 600 400]);
        decision=1; %1-hard decision, 0-soft decision
        [est_info_bits_MFGS,det_iters_MFGS,~] = Matched_Filter_GS_detector(N,M,M_mod,sigma,data_grid,Y,H_t_f,n_ite_MRC,omega,Tn_block_matrix,Gn_block_matrix,zn_block_vector,r,Wn,decision);
        [est_info_bits_1tap,~] = TF_single_tap_equalizer(N,M,M_mod,sigma,data_grid,Y,H_t_f,Wn);
        [est_info_bits_LMMSE,~] = Block_LMMSE_detector(N,M,M_mod,sigma,data_grid,Gn_block_matrix,r,Wn);
       
