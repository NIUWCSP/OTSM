function [Rx] = PlutoSet(txdata)   

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
    
% System Object Configuration
    s = iio_sys_obj_matlab; % MATLAB libiio Constructor
    s.ip_address = ip;
    s.dev_name = 'ad9361';
    s.in_ch_no = 2;
    s.out_ch_no = 2;
    s.in_ch_size = length(txdata);
    s.out_ch_size = length(txdata) * 10;
    
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

    input{1} = real(txdata);
    input{2} = imag(txdata);
    output = stepImpl(s, input);%調用pluto的通道資料

    I = output{1};
    Q = output{2};
    Rx = I+1i*Q;
end
