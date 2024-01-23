%% **********************���d�ҶȾA�Ω����ۦ��۵o�ϥ�-������********************************
clearvars -except times;close all;warning off; %�w�]����
set(0,'defaultfigurecolor','w'); 
%�[�Jpath
addpath ..\..\library 
addpath ..\..\library\matlab 
addpath ..\..\code\matlab\OFDM

%�R��.mat
if(0)
    Delete_mat;
end

%�]�w�������榸��
set_looptimes = 10;

% ���J���
if exist('ScattorData.mat','file')
    load('ScattorData.mat');
else
    AllRxDataSymbEqAverage = zeros(0,0); %�x�sScattor���
end
if exist('BERData.mat','file')
    load('BERData.mat');
else
    AllBERData = zeros(0,0); %�x�sBER���
end

%�]�wpluto IP
ip = '192.168.2.1';

%�]�w�P�i�JTX�禡
upsample=4; %�L���˨�4���A�Ʀ��٭���������i�H�����u
txdata = Transmitter(upsample);
txdata = round(txdata .* 2^15); %[(32768*4096)/1024]  -2/2



%% Transmit and Receive using MATLAB libiio �걵pluto

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



%% Initializing simulation error count variables

N_fram = 1000;

est_info_bits_MFGS=zeros(N_bits_perfram,1);
est_info_bits_1tap=zeros(N_bits_perfram,1);
est_info_bits_LMMSE=zeros(N_bits_perfram,1);


err_ber_MFGS = zeros(1,set_looptimes);%bit error rate
err_ber_1tap = zeros(1,set_looptimes);
err_ber_LMMSE = zeros(1,set_looptimes);

avg_ber_MFGS=zeros(1,set_looptimes);
avg_ber_1tap=zeros(1,set_looptimes);
avg_ber_LMMSE=zeros(1,set_looptimes);

det_iters_MFGS=0;
no_of_detetor_iterations_MFGS= zeros(1,set_looptimes);
avg_no_of_iterations_MFGS=zeros(1,set_looptimes);

%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix
current_frame_number=zeros(1,set_looptimes);

%% PLOT TX for Evan_debug �e�XTX���ɰ�ϻP�W�й�
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
    
spectrum1 = dsp.SpectrumAnalyzer('SampleRate',      samp_rate, ...
                                'SpectrumType',    'Power density', ...
                                'SpectralAverages', 10, ...
                                'YLimits',         [-70 60], ...
                                'Title',           SpectrumTitleStr, ...
                                'YLabel',          'Power spectral density', ...
                                'Position',        [500 300 400 400]);

% Show power spectral density of captured burst
step(spectrum1,txdata);
release(spectrum1);


global NoFoundDataTimes;
NoFoundDataTimes = 0; %�����data������(�S���q�L�P�B)
global AllFoundDataTimes;
AllFoundDataTimes= 0; %���ձ���data������(�]�O������榸��)
AllRxDataSymbEq = zeros(0,0); %�Ȯ��x�sscattor���
AllBERDataCol = zeros(0,0); %�Ȯ��x�sBER���

for loop_times = 1:set_looptimes %�T�O���ƺʴ�
%% PLOT RX �e�XRX����
for i=i:4 %�ѩ�PLUTO-USB�ƾڶq����~�]��RX�ϥΦ�FOR-LOOP����TX�ƾڶi�J by Evan 2019-04-16
    fprintf('Transmitting Data Block %i ...\n',i);
    input{1} = real(txdata);
    input{2} = imag(txdata);
    output = stepImpl(s, input);
    fprintf('Data Block %i Received...\n',i);
end
    I = output{1};
    Q = output{2};
    Rx = I+1i*Q;
    figure(1); clf;
    set(gcf,'name','������-RX���I/Q�������A'); % EVAN for debug OK
    subplot(121);
    plot(I);
    hold on;
    plot(Q);
    subplot(122);
    pwelch(Rx, [],[],[], 40e6, 'centered', 'psd');
    % 20230301�s�W�NPSD���|�_��
    hold on; 
    pwelch(txdata, [],[],[], 40e6, 'centered', 'psd');
    legend('Rx', 'Tx')
    %% PLOT RX
    R6x = Rx(:,1);
    global RxDataSymbEq;
    Receiver(Rx(1:4:end));
    % �x�sscattor���(�@�����108*18�ҥH�C18�C�@��)
    AllRxDataSymbEq = [AllRxDataSymbEq RxDataSymbEq];
    

    global TxDataBits;
    global RxDataBits;
    global BERData;
    AllFoundDataTimes = AllFoundDataTimes + 1;
    BER ( TxDataBits , RxDataBits , NoFoundDataTimes , AllFoundDataTimes);
    AllBERDataCol = [AllBERDataCol BERData];
    
    %pause(0.1);
    %break;
end

%% ����
fprintf('Transmission and reception finished\n');

% Read the RSSI attributes of both channels
rssi1 = output{s.getOutChannel('RX1_RSSI')};
% rssi2 = output{s.getOutChannel('RX2_RSSI')};

s.releaseImpl();

% �����Y�@�T�w�Z����BER��P�y��
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

% �N��ƾɥX��.mat��
save('ScattorData.mat', 'AllRxDataSymbEqAverage');
save('BERData.mat', 'AllBERData');
