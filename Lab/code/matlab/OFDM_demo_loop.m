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

%設定接收執行次數
global set_looptimes 
set_looptimes= 5;

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

%% Initializing simulation error count variables

N_fram = 1000;

err_ber_MFGS = zeros(1,set_looptimes);%bit error rate
err_ber_1tap = zeros(1,set_looptimes);
err_ber_LMMSE = zeros(1,set_looptimes);

avg_ber_MFGS=zeros(1,set_looptimes);
avg_ber_1tap=zeros(1,set_looptimes);
avg_ber_LMMSE=zeros(1,set_looptimes);

det_iters_MFGS=0;
no_of_detetor_iterations_MFGS= zeros(set_looptimes,1); %no_of_detetor_iterations_MFGS= zeros(1,set_looptimes);
avg_no_of_iterations_MFGS=zeros(1,set_looptimes);

%% Transmit and Receive using MATLAB libiio 串接pluto
        
% System Object Configuration
 s = iio_sys_obj_matlab; % MATLAB libiio Constructor
 s.ip_address = ip;
 s.dev_name = 'ad9361';
 s.in_ch_no = 2;
 s.out_ch_no = 2;
 s.in_ch_size = N*M;
 s.out_ch_size = N*M.*10;
        
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
         
        
        
        %% PLOT TX for Evan_debug 畫出TX的時域圖與頻譜圖
        
        TimeScopeTitleStr = 'OFDM-TX-Baseband I/Q Signal';
        SpectrumTitleStr = 'OFDM-TX-Baseband Signal Spectrum';
            
        samp_rate = 40e6;    
            scope1 = dsp.TimeScope('SampleRate',     samp_rate, ...
                                  'Title',           TimeScopeTitleStr, ...
                                  'TimeSpan',        1/samp_rate*(length(txdata)+100), ...
                                  'YLimits',         [-5 5], ...
                                  'ShowLegend',      true, ...
                                  'ShowGrid',        true, ...
                                  'Position',        [0 300 400 400]);
                                  
            step(scope1,txdata);
            release(scope1);
        
        
        global NoFoundDataTimes;
        NoFoundDataTimes = 0; %未找到data的次數(沒有通過同步)
        global AllFoundDataTimes;
        AllFoundDataTimes= 0; %嘗試接收data的次數(也是環圈執行次數)
        AllRxDataSymbEq = zeros(0,0); %暫時儲存scattor資料
        AllBERDataCol = zeros(0,0); %暫時儲存BER資料
        
%% Initializing simulation error count variables
N_fram = 1000;
for iesn0 = 1:set_looptimes
    for ifram = 1:N_fram 
        current_frame_number=zeros(1,iesn0);
        current_frame_number(iesn0)=ifram;

        %% PLOT RX 畫出RX的圖
        for i=i:4 %由於PLUTO-USB數據量受限~因此RX使用此FOR-LOOP等待TX數據進入 by Evan 2019-04-16
            fprintf('Transmitting Data Block %i ...\n',i);
            input{1} = reshape(real(txdata),[],1);
            input{2} = reshape(imag(txdata),[],1);
            output = stepImpl(s, input);%調用pluto的通道資料
            fprintf('Data Block %i Received...\n',i);
        end
            I = output{1};
            Q = output{2};
            Rx = I+1i*Q;
            figure(1); clf;
            set(gcf,'name','立鎂科技-RX實際I/Q接收狀態'); % EVAN for debug OK
            subplot(121);
            plot(I);
            hold on;
            plot(Q);
            subplot(122);
            pwelch(Rx, [],[],[], 40e6, 'centered', 'psd');
            % 20230301新增將PSD圖疊起來
            hold on; 
            pwelch(txdata, [],[],[], 40e6, 'centered', 'psd');
            legend('Rx', 'Tx')
            %% PLOT RX
            R6x = Rx(:,1);
            global RxDataSymbEq;
            [RxDataBits,est_info_bits_MFGS,det_iters_MFGS,est_info_bits_1tap,est_info_bits_LMMSE] = Receiver(Rx(1:4:end));
            % 儲存scattor資料(一筆資料108*18所以每18列一筆)
            AllRxDataSymbEq = [AllRxDataSymbEq RxDataSymbEq];
            
        
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
        display(set_looptimes,'set_looptimes');
        display(avg_ber_MFGS,'Average BER - Matched Filtered Gauss Seidel');
        display(avg_ber_1tap,'Average BER - single tap equalizer');
        display(avg_ber_LMMSE,'Average BER - LMMSE equalizer');
        display(avg_no_of_iterations_MFGS,'Average number of iterations for the MFGS detector');
        disp('####################################################################')
      end
end


%% 結束
fprintf('Transmission and reception finished\n');

% Read the RSSI attributes of both channels
rssi1 = output{s.getOutChannel('RX1_RSSI')};
% rssi2 = output{s.getOutChannel('RX2_RSSI')};

s.releaseImpl();

% 平均某一固定距離的BER跟星座圖
AllBERData = [AllBERData;AllBERDataCol];
for i=1:108
    for j=1:18
        for k=1:set_looptimes
            AllRxDataSymbEq(i,j) = AllRxDataSymbEq(i,j)+AllRxDataSymbEq(i,j+18*(k-1));
            
        end
    end
end
AverageData = zeros(108,18);
AverageData = AllRxDataSymbEq(1:108,1:18)./set_looptimes;
AllRxDataSymbEqAverage = [AllRxDataSymbEqAverage AverageData];

% 將資料導出到.mat中
save('ScattorData.mat', 'AllRxDataSymbEqAverage');
save('BERData.mat', 'AllBERData');
