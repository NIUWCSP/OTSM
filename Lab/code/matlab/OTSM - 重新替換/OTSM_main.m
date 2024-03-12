%% **********************此範例僅適用於單機自收自發使用-立鎂科技********************************
upsample=4; %過取樣取4倍，數位還原類比後比較可以不失真
[AllRxDataSymbEqAverage, AllBERData,ip, upsample] = pluto();
%txdata = Transmitter(upsample);
%txdata = round(txdata .* 2^15);

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
%delay_spread = M/(8/3);%40*64是資料部分 剩下是Pilot跟Sync
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

N_fram = 100;

est_info_bits_MFGS=zeros(N_bits_perfram,1);
est_info_bits_1tap=zeros(N_bits_perfram,1);
est_info_bits_LMMSE=zeros(N_bits_perfram,1);


err_ber_MFGS = zeros(1,length(SNR_dB));%bit error rate
err_ber_1tap = zeros(1,length(SNR_dB));
err_ber_LMMSE = zeros(1,length(SNR_dB));

avg_ber_MFGS=zeros(1,length(SNR_dB));
avg_ber_1tap=zeros(1,length(SNR_dB));%set_looptimes
avg_ber_LMMSE=zeros(1,length(SNR_dB));

det_iters_MFGS=0;
no_of_detetor_iterations_MFGS= zeros(length(SNR_dB),1);
avg_no_of_iterations_MFGS=zeros(1,length(SNR_dB)); 

% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix
current_frame_number=zeros(1,length(SNR_dB));
%% Transmitter
[tx_signal2,TxDataBits] = Transmitter(upsample,N,M,M_mod,M_bits,data_grid,N_syms_perfram,Wn);
%% Transmit and Receive using MATLAB libiio 串接pluto

[input, output,s] = configureAD9361(ip, txdata); % System Object Configuration

for iesn0 = 1:length(SNR_dB)  %iesn0=loop_times 
    for ifram = 1:N_fram
        current_frame_number(iesn0)=ifram;
        %% random input bits generation%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);%trans_info_bit =TxDataBits
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');%data=1*3840        
        %data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod, 0,'gray','bit');  data=2*3840      
        X = Generate_2D_data_grid(N,M,data,data_grid);
        
        
        %% OTSM modulation%%%%
        X_tilda=X*Wn;               %equation (6) in [R1]
        s = reshape(X_tilda,N*M,1); %equation (7) in [R1]
        
        
        %% OTFS channel generation%%%%
        
        %         %% test channel
        %         taps=4;
        %         chan_coef=1/sqrt(2)*(randn(1,taps)+1i.*randn(1,taps));
        %         delay_taps=[0,1,2,3];  %% maximum value should be less than the channel delay spread
        %         Doppler_taps=[0,1,2,3];
        
        % 3GPP channel model
        max_speed=500;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        
        
        
        
        %% channel output%%%%%
        [G,gs]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);
        
        r=zeros(N*M,1);
        noise= sqrt(sigma_2(iesn0)/2)*(randn(size(s)) + 1i*randn(size(s)));
        l_max=max(delay_taps);
        for q=0:N*M-1
            for l=0:l_max
                if(q>=l)
                    r(q+1)=r(q+1)+gs(l+1,q+1)*s(q-l+1);  %equation (24) in [R1]
                end
            end
        end
        r=r+noise;
        
        
        %% OTSM demodulation%%%%
        Y_tilda=reshape(r,M,N);     %equation (11) in [R1]
        Y = Y_tilda*Wn;             %equation (12) in [R1]
        
        
        %% test: the received time domain signal can be generated element by element (using gs) or in the matrix form (using r=G.s).
        %         r_test=G*s+noise;                            %equation (31) in [R1]
        %         test_delay_time_matrix_error=norm(r_test-r)
        %% test: the received delay-seq signal can be in the matrix form (using y=H.x).
        %         [H,P]= Gen_delay_sequency_channel_matrix(N,M,G,Wn);
        %         noise_DS=kron(eye(M),Wn)*P'*noise;            %equation (34)in [R1]
        %         x_vec=reshape(X,N*M,1);
        %         y_vec=reshape(Y,N*M,1);
        %         y_vec_test=H*x_vec+noise_DS;                %equation (33) in   [R1]
        %         text_delay_seq_matrix_error=norm(y_vec_test-y_vec)
        
        
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
        
        
        %% errors count%%%%%
        errors_MFGS = sum(xor(est_info_bits_MFGS,trans_info_bit));
        errors_1tap = sum(xor(est_info_bits_1tap,trans_info_bit));
        errors_LMMSE = sum(xor(est_info_bits_LMMSE,trans_info_bit));
        
        
        err_ber_MFGS(1,iesn0) = err_ber_MFGS(1,iesn0) + errors_MFGS;
        err_ber_1tap(1,iesn0) = err_ber_1tap(1,iesn0) + errors_1tap;
        err_ber_LMMSE(1,iesn0) = err_ber_LMMSE(1,iesn0) + errors_LMMSE;
        
        no_of_detetor_iterations_MFGS(iesn0)=no_of_detetor_iterations_MFGS(iesn0)+det_iters_MFGS;
        
        
        
        %%  Error count
        
        avg_no_of_iterations_MFGS(iesn0)=no_of_detetor_iterations_MFGS(iesn0)/ifram;
        avg_ber_MFGS(1,iesn0)=err_ber_MFGS(1,iesn0).'/length(trans_info_bit)/ifram;
        avg_ber_1tap(1,iesn0)=err_ber_1tap(1,iesn0).'/length(trans_info_bit)/ifram;
        avg_ber_LMMSE(1,iesn0)=err_ber_LMMSE(1,iesn0).'/length(trans_info_bit)/ifram;
        
        
        
        %%         DISP error performance details
        clc
        disp('####################################################################')
        fprintf('OTSM-(N,M,QAM size)');disp([N,M,M_mod]);
        display(current_frame_number,'Number of frames');
        display(SNR_dB,'SNR (dB)');
        display(avg_ber_MFGS,'Average BER - Matched Filtered Gauss Seidel');
        display(avg_ber_1tap,'Average BER - single tap equalizer');
        display(avg_ber_LMMSE,'Average BER - LMMSE equalizer');
        display(avg_no_of_iterations_MFGS,'Average number of iterations for the MFGS detector');
        disp('####################################################################')
        
    end
    
end
figure(1)
semilogy(SNR_dB,avg_ber_MFGS,'-o','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_1tap,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_LMMSE,'-s','LineWidth',2,'MarkerSize',8)
legend('MFGS','single tap','LMMSE')
