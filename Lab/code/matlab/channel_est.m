function [ChanEst] = channel_est()

XCorrPilot = Y_OTSM_Pilot(:,1)' * Y_OTSM_Pilot(:,2);%導頻符號的交叉相關
EpsEst = 1/(2*pi) * atan(imag(XCorrPilot)/real(XCorrPilot));%頻率偏移的估計

% Estimate carrier freqnecy offset
RxSigalRadioFrameCmpCFO = RxSignalRadioFrame .* ...
    exp(-1j*2*pi*EpsEst/NumFFT * (0:length(RxSignalRadioFrame)-1)');%接收訊號進行CFO校正
RxSignalRadioGridCFO = reshape(RxSigalRadioFrameCmpCFO,N,M);

% Reobtain pilot data
PilotOtsmSymb = RxSignalRadioFrame(M_data+sqrt(PilotSymb)+1:M_data+sqrt(PilotSymb)*2,1:sqrt(PilotSymb));%將所選的一段導頻資料重新組織成一個矩陣，其中每列有兩個元素
PilotOtsmSymb = RxSignalRadioGridCFO(M_data+sqrt(PilotSymb)+1:M_data+sqrt(PilotSymb)*2,1:sqrt(PilotSymb));%將所選的一段導頻資料重新組織成一個矩陣，其中每列有兩個元素

% Data OFDM symbol
DataOtsmSymb = RxSignalRadioGridCFO(1:M_data,1:M);


% Channel estimation and equalization
PilotBits = GetPilotBits();
Tx_PilotSymb = qammod(reshape(PilotBits,M_bits,[]),M_mod,'gray','InputType','bit');%產生Tx的 BPSK符號序列
ChanEst = reshape(PilotOtsmSymb,1,[]) ./ Tx_PilotSymb;%通道估计
global RxDataSymbEq;
RxDataSymbEq = DataOtsmSymb ./ repmat(ChanEst, M_data,1);
subplot(232);plot(10*log10(abs(ChanEst).^2)-min(10*log10(abs(ChanEst).^2)));title('channel estimation');%繪製通道估計的幅度譜
subplot(233);plot(DataOtsmSymb(:),'*');axis equal;title('scatter before equalization');axis square;
subplot(234);plot(RxDataSymbEq(:).*exp(-1i*pi/4),'.');axis([-1.5,1.5,-1.5,1.5]);title('scatter after equalization'); axis square;%*exp(-1i*pi/4) 的作用是進行相位調整

end