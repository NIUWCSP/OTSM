function tx_signal2 = Transmitter(upsample,N,M,M_mod)
%NumFFT = 64;%V3 FFT轉換的點數
NumSyncPreamble = 32;%V3 同步的前綴，Preamble：防干擾+同步+通道估測(已知的頻域資料)
NumCP = 16;%V3 CP：循環前綴，CP：避免ISI(多路徑干擾)(未知的時域訊號)

%% OTFS parameters%%%%%%%%%%

M_bits = log2(M_mod);

%% delay-Doppler grid symbol placement
% max delay spread in the channel
delay_spread = M/(8/3);%40*64是資料部分 剩下是Pilot跟Sync
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-delay_spread;
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));

%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix

%% Transmitter
% Generate pilot symbols
PilotBits = GetPilotBits();%Preamble的data

% Generate synchronization symbols
SyncBits = GetSyncBits();
QamSyncBits=reshape(qammod(reshape(SyncBits,M_bits,size(SyncBits,2)/2), M_mod,'gray','InputType','bit'),[],1);
QamSync_tilda = CompeteISI(QamSyncBits,M,0); %128*1

% Generate data symbols
global TxDataBits;
TxDataBits = GetTxDataBits();%TX的data
TxData=qammod(reshape(TxDataBits,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');%data=1*2560
Tx = Generate_2D_data_grid(N,M,TxData,data_grid);
Tx_Symb=Tx_addPilot(Tx,PilotBits,N,M_mod); 

%% OTSM modulation%%%%
Tx_tilda=Tx_Symb*Wn;              %equation (6) in [R1]   %Tx=X
%TxDataOtsmSymb=OtsmSignalModulation(Tx_tilda, NumFFT, NumCP);
tx_Data_signal=CompeteISI(reshape(Tx_tilda,[],1),M,NumCP);  %equation (7) in [R1]
TxSignal = [ ...
    QamSync_tilda(1:NumSyncPreamble);
    QamSync_tilda(1:NumSyncPreamble);
    QamSync_tilda;
    tx_Data_signal];
flt1=rcosine(1,upsample,'fir/sqrt',0.05,64);%pulse shaper 
tx_signal2=rcosflt(TxSignal,1,upsample, 'filter', flt1); %4240(TxSignal)*4(upsample)+(513(flt1)-1)：因為捲積所以要-1
%tx_signal2=TxSignal;

%%陳昱升想測試畫圖%%

figure(1);
set(gcf,'name','Transmitter內的資料實際點數','Position', [20 400 500 500]);
subplot(221);
plot(TxDataBits,'.');title('Orignal Txdatabits');axis([-10,N_syms_perfram*M_bits+10,-0.5,1.5]);
subplot(222);
plot(TxData,'*');title('After Qammod Txdatabits');axis([-1.5,1.5,-1.5,1.5]);
subplot(223);
plot(Tx_tilda,'*');title('Then addPilot&WHT');
subplot(224);
plot(tx_signal2,'.');title('Final TxSignal');
