function [RxDataSymbEq] = channel_est(N,M,M_mod,NumFFT,RxSignalRadioFrame,Y_OTSM_Pilot)

%%基本參數設置
% max delay spread in the channel
delay_spread = M/(8/3);%40*64是資料部分 剩下是Pilot跟Sync
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-delay_spread;%64-24=40
M_bits = log2(M_mod);
NumDataN=N-delay_spread;
PilotSymb=size(GetPilotBits,2)/2;%PilotBits為128個Bits QAM後會除2

%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix

% Estimate carrier frequency offset
XCorrPilot = Y_OTSM_Pilot(:,1)' * Y_OTSM_Pilot(:,2);%導頻符號的交叉相關
EpsEst = 1/(2*pi) * atan(imag(XCorrPilot)/real(XCorrPilot));%頻率偏移的估計

% Estimate carrier freqnecy offset
RxSigalRadioFrameCmpCFO = RxSignalRadioFrame .* ...
    exp(-1j*2*pi*EpsEst/NumFFT * (0:length(RxSignalRadioFrame)-1)');%接收訊號進行CFO校正
RxSignalRadioGridCFO = reshape(RxSigalRadioFrameCmpCFO,N,M);

% Reobtain pilot data
PilotOtsmSymb = RxSignalRadioGridCFO(M_data+sqrt(PilotSymb)+1:M_data+sqrt(PilotSymb)*2,1:sqrt(PilotSymb));%將所選的一段導頻資料重新組織成一個矩陣，其中每列有兩個元素

% Data OFDM symbol
DataOtsmSymb = RxSignalRadioGridCFO(1:M_data,1:M);


% Channel estimation and equalization
PilotBits = GetPilotBits();
QamPilotSymb = qammod(reshape(PilotBits,M_bits,[]),M_mod,'gray','InputType','bit');
QamPilotSymbGrid= reshape(QamPilotSymb,sqrt(size(QamPilotSymb,2)),[]);
WnPilotSymb=zeros(N,M);
WnPilotSymb(NumDataN+sqrt(size(QamPilotSymb,2))+1:NumDataN+sqrt(size(QamPilotSymb,2))*2,1:sqrt(size(QamPilotSymb,2))) = QamPilotSymbGrid;
WnPilotSymb=WnPilotSymb*Wn;
Tx_PilotSymb=reshape(WnPilotSymb(NumDataN+sqrt(size(QamPilotSymb,2))+1:NumDataN+sqrt(size(QamPilotSymb,2))*2,1:sqrt(size(QamPilotSymb,2))),1,[]);

ChanEst = reshape(PilotOtsmSymb,1,[]) ./ Tx_PilotSymb;%通道估计
global RxDataSymbEq;
RxDataSymbEq = DataOtsmSymb ./ repmat(ChanEst, N,1);
subplot(232);plot(10*log10(abs(ChanEst).^2)-min(10*log10(abs(ChanEst).^2)));title('channel estimation');%繪製通道估計的幅度譜
subplot(233);plot(DataOtsmSymb(:),'*');axis equal;title('scatter before equalization');axis square;
subplot(234);plot(RxDataSymbEq(:).*exp(-1i*pi/4),'.');axis([-1.5,1.5,-1.5,1.5]);title('scatter after equalization'); axis square;%*exp(-1i*pi/4) 的作用是進行相位調整

end