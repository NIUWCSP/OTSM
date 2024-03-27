function [RxDataSymbEq,RxSigalRadioFrameCmpCFO,G] = channel_est(N,M,M_mod,RxSignalRadioFrame,Y_OTSM_Pilot,l_max)

%%基本參數設置
% max delay spread in the channel
delay_spread = M/(8/3);%40*64是資料部分 剩下是Pilot跟Sync
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-delay_spread;%64-24=40
M_bits = log2(M_mod);
DelayPilotSymb=sqrt(size(GetPilotBits,2)/2);%PilotBits為128個Bits QAM後會除2
RxSignalRadioGrid=reshape(RxSignalRadioFrame,N,M);

%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix


               
%% Gen_time_domain_channel

% Estimate carrier frequency offset
XCorrPilot = Y_OTSM_Pilot(:,1)' * Y_OTSM_Pilot(:,2);%導頻符號的交叉相關
EpsEst = 1/(2*pi) * atan(imag(XCorrPilot)/real(XCorrPilot));%頻率偏移的估計

 %測試用CFO
y_mp1=zeros(size(Y_OTSM_Pilot,1),l_max+1); %64*1
y_mp2=zeros(size(Y_OTSM_Pilot,1),l_max+1); %64*1
N_p2=sqrt(size(Y_OTSM_Pilot,1));

  for l=0:l_max
      for n=1:DelayPilotSymb %1~8
          for mp1=M_data+1:M_data+DelayPilotSymb %41~48 刪gs(l+1,mp1+l+n*M+1).*
                    y_mp1((mp1-M_data)+(n-1)*8,l+1)= exp(1j*2*pi*EpsEst/M*(mp1+l+n*M+1))*RxSignalRadioGrid(mp1,n);
          end
          for mp2=M_data+N_p2+1:M_data+DelayPilotSymb*2 %49~56 刪exp(1j*pi/2)*gs(l+1,mp2+l+n*M+1).*
                    y_mp2((mp2-M_data-N_p2)+(n-1)*8,l+1)= exp(1j*2*pi*EpsEst/M*(mp2+l+n*M+1))*RxSignalRadioGrid(mp2,n);
          end
      end
  end %%equation (17) in [R3]
    tilda_CFO=0; %初始化
  for l=0:l_max %把tilda_CFO的l_max*length(delay_taps)→*1
    XCorrmp=y_mp1(:,l+1) * (exp(1j*pi/2)*y_mp2(:,l+1)');
    tilda_CFO = tilda_CFO+(-M/(2*pi*1)* atan(imag(XCorrmp(:,l+1))/real(XCorrmp(:,l+1)))); %%equation (18) in [R3]
  end


% Estimate carrier freqnecy offset
 RxSigalRadioFrameCmpCFO = RxSignalRadioGrid.*exp(-1j*2*pi*tilda_CFO/M) ; %%equation (19) in [R3]
% RxSigalRadioFrameCmpCFO = RxSignalRadioFrame .* ...
%     exp(-1j*2*pi*EpsEst/M * (0:length(RxSignalRadioFrame)-1)');%接收訊號進行CFO校正
% RxSigalRadioFrameCmpCFO=reshape(RxSigalRadioFrameCmpCFO,N,M);
PilotOtsmSymb = RxSigalRadioFrameCmpCFO(M_data+1:M_data+DelayPilotSymb*1,1:DelayPilotSymb);%將所選的一段導頻資料重新組織成一個矩陣，其中每列有兩個元素


% Data OFDM symbol
DataOtsmSymb = RxSigalRadioFrameCmpCFO(1:M_data,1:M);

% Channel estimation and equalization
PilotBits = GetPilotBits();
QamPilotSymb = qammod(reshape(PilotBits,M_bits,[]),M_mod,'gray','InputType','bit');
QamPilotSymbGrid= reshape(QamPilotSymb,DelayPilotSymb,[]);
WnPilotSymb=zeros(N,M);
WnPilotSymb(M_data+1:M_data+DelayPilotSymb,1:DelayPilotSymb) = QamPilotSymbGrid;
WnPilotSymb=WnPilotSymb*Wn;
Tx_PilotSymb=WnPilotSymb(M_data+1:M_data+DelayPilotSymb,1:DelayPilotSymb);

ChanEst = PilotOtsmSymb ./ Tx_PilotSymb;%通道估计

%%試求gs、G
G=zeros(N*M,N*M);
gs_Grid=repmat(ChanEst, N/size(ChanEst,1),M/size(ChanEst,2));
gs=reshape(gs_Grid,l_max+1,[]);

for q=0:N*M-1
    for l=0:l_max
        if(q>=l)
            G(q+1,q+1-l)=gs(l+1,q+1);
        end
    end
end

RxDataSymbEq = RxSigalRadioFrameCmpCFO./gs_Grid;

%%偵錯
global NoFoundDataTimes;
if (isnan(atan(imag(XCorrmp)/real(XCorrmp))))
    NoFoundDataTimes = NoFoundDataTimes+1;
    RxDataSymbEq(1:N,1:M)=RxSignalRadioGrid(1:N,1:M);
    return;
end

% RxDataSymbEq = DataOtsmSymb ./ repmat(ChanEst, M_data/size(ChanEst,1),M/size(ChanEst,2));
subplot(232);plot(10*log10(abs(reshape(ChanEst,[],1)).^2)-min(10*log10(abs(reshape(ChanEst,[],1)).^2)));title('channel estimation');%繪製通道估計的幅度譜
subplot(233);plot(DataOtsmSymb(:),'*');axis equal;title('scatter before equalization');axis square;
subplot(234);plot(RxDataSymbEq(:),'.');axis([-1.5,1.5,-1.5,1.5]);title('scatter after equalization'); axis square;%*exp(-1i*pi/4) 的作用是進行相位調整
subplot(235);plot([RxDataSymbEq;zeros(N-M_data,M)]*Wn,'.');axis([-1.5,1.5,-1.5,1.5]);title('After WHT'); axis square;
text(2,0.6,['NoFoundData: ',num2str(NoFoundDataTimes),' 次']);

end