function [RxDataBits,est_info_bits_MFGS,det_iters_MFGS,est_info_bits_1tap,est_info_bits_LMMSE] = Receiver(RxSignal)

%%廣域變數宣告
global iesn0
global ifram

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
delay_spread = M/16;
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


%% Initializing simulation error count variables
est_info_bits_MFGS=zeros(N_bits_perfram,1);
est_info_bits_1tap=zeros(N_bits_perfram,1);
est_info_bits_LMMSE=zeros(N_bits_perfram,1);

%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix


         %% OTFS channel generation%%%%
         % 3GPP channel model
         max_speed=500;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        
         %% channel output%%%%% 
        [G,gs]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);
        
         r=zeros(N*M,1);
         noise= sqrt(sigma_2(iesn0)/2)*(randn(size(RxSignal)) + 1i*randn(size(RxSignal)));
         l_max=max(delay_taps);
         for q=0:N*M-1
           for l=0:l_max
             if(q>=l)
               r(q+1)=r(q+1)+gs(l+1,q+1)*RxSignal(q-l+1);  %equation (24) in [R1]
             end
           end
         end
         r=r+noise;

        %% OTSM demodulation%%%%
                Y_tilda=reshape(r,M,N);     %equation (11) in [R1]
                Y = Y_tilda*Wn;             %equation (12) in [R1]

        %% OTSM Reget Pilot%%%%
                Y_tilda_Pilot=reshape(Y_tilda(size(Y_tilda,1),:),[],1);
                PilotBits = GetPilotBits();
        
        %% Generate the block-wise channel matrices in the delay-time and the time-frequency domain
        [Gn_block_matrix,Tn_block_matrix,zn_block_vector,H_t_f]=Generate_Matched_Filter_GS_matrices(N,M,G,r);
        
         %% GS SOR Iterative detection
        
        n_ite_MRC=50; % maximum number of detector iterations
        omega=1; %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        if(M_mod==64)
            omega=0.25;
        end
        decision=1; %1-hard decision, 0-soft decision
        [est_info_bits_MFGS,det_iters_MFGS,data_MFGS] = Matched_Filter_GS_detector(N,M,M_mod,sigma_2(iesn0),data_grid,Y,H_t_f,n_ite_MRC,omega,Tn_block_matrix,Gn_block_matrix,zn_block_vector,r,Wn,decision);
        [est_info_bits_1tap,data_1tap] = TF_single_tap_equalizer(N,M,M_mod,sigma_2(iesn0),data_grid,Y,H_t_f,Wn);
        [est_info_bits_LMMSE,data_LMMSE] = Block_LMMSE_detector(N,M,M_mod,sigma_2(iesn0),data_grid,Gn_block_matrix,r,Wn);
        RxDataBits=0;
       
