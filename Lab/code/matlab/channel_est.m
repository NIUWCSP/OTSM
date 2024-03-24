function [RxDataSymbEq] = channel_est(N,M,M_mod,RxSignalRadioFrame,Y_OTSM_Pilot,l_max,gs,delay_taps)

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

%測試用CFO
RxSignalRadioGrid=reshape(RxSignalRadioFrame,N,M);
y_mp1=zeros(size(Y_OTSM_Pilot,1),1); %64*1
y_mp2=zeros(size(Y_OTSM_Pilot,1),1); %64*1

  for l=0:l_max
      for n=1:sqrt(size(Y_OTSM_Pilot,1)) %1~8
          for mp1=M_data+1:M_data+sqrt(size(Y_OTSM_Pilot,1)) %41~48
                    y_mp1((mp1-M_data)+(n-1)*8,1)=gs(l+1,mp1+l+n*M+1).* exp(1j*2*pi*EpsEst/M*(mp1+l+n*M+1))*RxSignalRadioGrid(mp1,n);
          end
          for mp2=M_data+sqrt(size(Y_OTSM_Pilot,1))+1:M_data+sqrt(size(Y_OTSM_Pilot,1))*2 %49~56
                    y_mp2((mp2-M_data-sqrt(size(Y_OTSM_Pilot,1)))+(n-1)*8,1)=gs(l+1,mp2+l+n*M+1).* exp(1j*2*pi*EpsEst/M*(mp2+l+n*M+1))*RxSignalRadioGrid(mp2,n);
          end
      end
  end %%equation (17) in [R3]
XCorrmp = y_mp1(:,1) * (exp(1j*pi/2)*y_mp2(:,1)');
tilda_CFO = -M/(2*pi*l_max*length(delay_taps))* atan(imag(XCorrmp)/real(XCorrmp)); %%equation (18) in [R3]


% Estimate carrier freqnecy offset
RxSigalRadioFrameCmpCFO = RxSignalRadioGrid.*exp(-1j*2*pi*tilda_CFO/M) ;

%RxSigalRadioFrameCmpCFO = RxSignalRadioFrame .* ...
    %exp(-1j*2*pi*EpsEst/M * (0:length(RxSignalRadioFrame)-1)');%接收訊號進行CFO校正
RxSignalRadioGridCFO = reshape(RxSigalRadioFrameCmpCFO,N,M);

% Reobtain pilot data
PilotOtsmSymb = (RxSignalRadioGridCFO(M_data+1:M_data+sqrt(PilotSymb)*1,1:sqrt(PilotSymb)));%將所選的一段導頻資料重新組織成一個矩陣，其中每列有兩個元素

% Data OFDM symbol
DataOtsmSymb = RxSignalRadioGridCFO(1:M_data,1:M);


% Channel estimation and equalization
PilotBits = GetPilotBits();
QamPilotSymb = qammod(reshape(PilotBits,M_bits,[]),M_mod,'gray','InputType','bit');
QamPilotSymbGrid= reshape(QamPilotSymb,sqrt(size(QamPilotSymb,2)),[]);
WnPilotSymb=zeros(N,M);
%WnPilotSymb(NumDataN+sqrt(size(QamPilotSymb,2))+1:NumDataN+sqrt(size(QamPilotSymb,2))*2,1:sqrt(size(QamPilotSymb,2))) = QamPilotSymbGrid;
WnPilotSymb(NumDataN+1:NumDataN+sqrt(size(QamPilotSymb,2)),1:sqrt(size(QamPilotSymb,2))) = QamPilotSymbGrid;
WnPilotSymb=WnPilotSymb*Wn;
%Tx_PilotSymb=WnPilotSymb(NumDataN+sqrt(size(QamPilotSymb,2))+1:NumDataN+sqrt(size(QamPilotSymb,2))*2,1:sqrt(size(QamPilotSymb,2)));
Tx_PilotSymb=WnPilotSymb(NumDataN+1:NumDataN+sqrt(size(QamPilotSymb,2)),1:sqrt(size(QamPilotSymb,2)));


ChanEst = PilotOtsmSymb ./ Tx_PilotSymb;%通道估计

RxDataSymbEq = DataOtsmSymb ./ repmat(ChanEst, M_data/size(ChanEst,1),M/size(ChanEst,2));
subplot(232);plot(10*log10(abs(reshape(ChanEst,[],1)).^2)-min(10*log10(abs(reshape(ChanEst,[],1)).^2)));title('channel estimation');%繪製通道估計的幅度譜
subplot(233);plot(DataOtsmSymb(:),'*');axis equal;title('scatter before equalization');axis square;
subplot(234);plot(RxDataSymbEq(:).*exp(-1i*pi/4),'.');axis([-1.5,1.5,-1.5,1.5]);title('scatter after equalization'); axis square;%*exp(-1i*pi/4) 的作用是進行相位調整
subplot(235);plot([RxDataSymbEq;zeros(N-M_data,M)]*Wn.*exp(-1i*pi/4),'.');axis([-1.5,1.5,-1.5,1.5]);title('After WHT'); axis square;

end