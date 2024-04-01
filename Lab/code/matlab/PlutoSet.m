function [Rx] = PlutoSet(txdata,sigma)   

%設定pluto IP
ip = '192.168.2.1';
    
% System Object Configuration
    s = iio_sys_obj_matlab; % MATLAB libiio Constructor
    s.ip_address = ip;
    s.dev_name = 'ad9361';
    s.in_ch_no = 2;
    s.out_ch_no = 2;
    s.in_ch_size = length(txdata);
    s.out_ch_size = length(txdata) * 4;
    
    s = s.setupImpl();
    
    input = cell(1, s.in_ch_no + length(s.iio_dev_cfg.cfg_ch));
    output = cell(1, s.out_ch_no + length(s.iio_dev_cfg.mon_ch));
    
    % Set the attributes of AD9361
    input{s.getInChannel('RX_LO_FREQ')} = 2400e6;
    input{s.getInChannel('RX_SAMPLING_FREQ')} = 40e6;
    input{s.getInChannel('RX_RF_BANDWIDTH')} = 20e6;
    input{s.getInChannel('RX1_GAIN_MODE')} = 'manual';%% slow_attack manual
    %input{s.getInChannel('TX1_GAIN')} = sqrt(sigma/2);
    input{s.getInChannel('RX1_GAIN')} = 1+sqrt(sigma/2);
    input{s.getInChannel('TX_LO_FREQ')} = 2400e6;
    input{s.getInChannel('TX_SAMPLING_FREQ')} = 40e6;
    input{s.getInChannel('TX_RF_BANDWIDTH')} = 20e6;
for i=0:4 %由於PLUTO-USB數據量受限~因此RX使用此FOR-LOOP等待TX數據進入 by Evan 2019-04-16
    input{1} = real(txdata);
    input{2} = imag(txdata);
    output = stepImpl(s, input);%調用pluto的通道資料
end
    I = output{1};
    Q = output{2};
    Rx = I+1i*Q;
    figure(2); clf;%clear figure
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

end
