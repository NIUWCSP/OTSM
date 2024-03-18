%% **********************此範例僅適用於單機自收自發使用-立鎂科技********************************
clearvars -except times;close all;warning off; %預設環境
set(0,'defaultfigurecolor','w'); 
%加入path
addpath ..\..\library 
addpath ..\..\library\matlab 
addpath ..\..\code\matlab\OFDM

%刪除.mat
if(0)
    Delete_mat;
end

% 載入資料
if exist('ScattorData.mat','file')
    load('ScattorData.mat');
else
    AllRxDataSymbEqAverage = zeros(0,0); %儲存Scattor資料
end
if exist('BERData.mat','file')
    load('BERData.mat');
else
    AllBERData = zeros(0,0); %儲存BER資料
end

%設定pluto IP
ip = '192.168.2.1';



%%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time

N = 64;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 10:2.5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;
%%%

%% Initializing simulation error count variables

err_ber_MFGS = zeros(1,length(SNR_dB));%bit error rate
err_ber_1tap = zeros(1,length(SNR_dB));
err_ber_LMMSE = zeros(1,length(SNR_dB));

avg_ber_MFGS=zeros(1,length(SNR_dB));
avg_ber_1tap=zeros(1,length(SNR_dB));
avg_ber_LMMSE=zeros(1,length(SNR_dB));

det_iters_MFGS=0;
no_of_detetor_iterations_MFGS= zeros(length(SNR_dB),1); %no_of_detetor_iterations_MFGS= zeros(1,set_looptimes);
avg_no_of_iterations_MFGS=zeros(1,length(SNR_dB));


%% Initializing simulation error count variables
N_fram = 10;
global iesn0
global ifram
for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram 
        current_frame_number=zeros(1,iesn0);
        current_frame_number(iesn0)=ifram;

        %設定與進入TX函式
        upsample=4; %過取樣取4倍，數位還原類比後比較可以不失真        
        txdata = Transmitter(upsample);
        txdata = round(txdata.*2^15);

%% PLOT RX 畫出RX的圖
Rx=PlutoSet(txdata);


            %% PLOT RX
            R6x = Rx(:,1);
            global RxDataSymbEq;
            [RxDataBits,est_info_bits_MFGS,det_iters_MFGS,est_info_bits_1tap,est_info_bits_LMMSE] = Receiver(txdata);
        
    %% errors count%%%%%
    global TxDataBits;
        errors_MFGS = sum(xor(est_info_bits_MFGS,TxDataBits));
        errors_1tap = sum(xor(est_info_bits_1tap,TxDataBits));
        errors_LMMSE = sum(xor(est_info_bits_LMMSE,TxDataBits));

        
        err_ber_MFGS(1,iesn0) = err_ber_MFGS(1,iesn0) + errors_MFGS;
        err_ber_1tap(1,iesn0) = err_ber_1tap(1,iesn0) + errors_1tap;
        err_ber_LMMSE(1,iesn0) = err_ber_LMMSE(1,iesn0) + errors_LMMSE;
        
        no_of_detetor_iterations_MFGS(iesn0)=no_of_detetor_iterations_MFGS(iesn0)+det_iters_MFGS;
        
        
        
        %%  Error count
        
        avg_no_of_iterations_MFGS(iesn0)=no_of_detetor_iterations_MFGS(iesn0)/ifram;
        avg_ber_MFGS(1,iesn0)=err_ber_MFGS(1,iesn0).'/length(TxDataBits)/ifram;
        avg_ber_1tap(1,iesn0)=err_ber_1tap(1,iesn0).'/length(TxDataBits)/ifram;
        avg_ber_LMMSE(1,iesn0)=err_ber_LMMSE(1,iesn0).'/length(TxDataBits)/ifram;
            
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


%% 結束
figure(3)
semilogy(SNR_dB,avg_ber_MFGS,'-o','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_1tap,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_LMMSE,'-square','LineWidth',2,'MarkerSize',8)
legend('MFGS','single tap','LMMSE')

