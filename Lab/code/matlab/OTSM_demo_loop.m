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

%設定與進入TX函式
upsample=4; %過取樣取4倍，數位還原類比後比較可以不失真        
txdata = Transmitter(upsample);

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

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 10:2.5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;

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

%% Transmit and Receive using MATLAB libiio 串接pluto
        
% System Object Configuration
 s = iio_sys_obj_matlab; % MATLAB libiio Constructor
 s.ip_address = ip;
 s.dev_name = 'ad9361';
 s.in_ch_no = 2;
 s.out_ch_no = 2;
 s.in_ch_size = length(txdata);
 s.out_ch_size = length(txdata).*10;
        
 s = s.setupImpl();
        
 input = cell(1, s.in_ch_no + length(s.iio_dev_cfg.cfg_ch));
 output = cell(1, s.out_ch_no + length(s.iio_dev_cfg.mon_ch));
        
 % Set the attributes of AD9361
  input{s.getInChannel('RX_LO_FREQ')} = 2400e6;
  input{s.getInChannel('RX_SAMPLING_FREQ')} = 40e6;
  input{s.getInChannel('RX_RF_BANDWIDTH')} = 20e6;
  input{s.getInChannel('RX1_GAIN_MODE')} = 'manual';%% slow_attack manual
  input{s.getInChannel('RX1_GAIN')} = 1;
  input{s.getInChannel('TX_LO_FREQ')} = 2400e6;
  input{s.getInChannel('TX_SAMPLING_FREQ')} = 40e6;
  input{s.getInChannel('TX_RF_BANDWIDTH')} = 20e6;
         

        
%% Initializing simulation error count variables
N_fram = 1000;
global iesn0
global ifram
for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram 
        current_frame_number=zeros(1,iesn0);
        current_frame_number(iesn0)=ifram;

%% PLOT RX 畫出RX的圖
for i=i:4 %由於PLUTO-USB數據量受限~因此RX使用此FOR-LOOP等待TX數據進入 by Evan 2019-04-16
    fprintf('Transmitting Data Block %i ...\n',i);
    input{1} = real(txdata);
    input{2} = imag(txdata);
    output = stepImpl(s, input);%調用pluto的通道資料
    fprintf('Data Block %i Received...\n',i);
end
    I = output{1};
    Q = output{2};
    Rx = I+1i*Q;
    figure(1); clf;%clear figure
    set(gcf,'name','立鎂科技-RX實際I/Q接收狀態'); % EVAN for debug OK %get current figure
    subplot(121);
    plot(I);
    hold on;
    plot(Q);
    subplot(122);
    pwelch(Rx, [],[],[], 40e6, 'centered', 'psd');
    % 20230301新增將PSD圖疊起來
    hold on; %'centered' 表示計算雙邊頻,'psd'表示頻譜類型
    pwelch(txdata, [],[],[], 40e6, 'centered', 'psd');
    legend('Rx', 'Tx')

            %% PLOT RX
            R6x = Rx(:,1);
            global RxDataSymbEq;
            [RxDataBits,est_info_bits_MFGS,det_iters_MFGS,est_info_bits_1tap,est_info_bits_LMMSE] = Receiver(Rx(1:10:end));
        
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
figure(2)
semilogy(SNR_dB,avg_ber_MFGS,'-o','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_1tap,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_LMMSE,'-s','LineWidth',2,'MarkerSize',8)
legend('MFGS','single tap','LMMSE')
