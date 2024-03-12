function [tx_signal2,TxDataBits] = Transmitter(upsample,N,M,M_mod,M_bits,data_grid,N_syms_perfram,Wn)
%NumFFT = 64;%V3 FFT轉換的點數
%NumSyncPreamble = 32;%V3 同步的前綴，Preamble：防干擾+同步+通道估測(已知的頻域資料)
NumCP = 16;%V3 CP：循環前綴，CP：避免ISI(多路徑干擾)(未知的時域訊號)

%% Transmitter
% Generate pilot symbols
PilotBits = GetPilotBits();%Preamble的data

% Generate synchronization symbols
SyncBits = GetSyncBits();%Preamble的data
%QamSyncBits=reshape(QamSyncBits,[],1);

% Generate data symbols
TxDataBits = randi([0,1],N_syms_perfram*M_bits,1);%TX的data
TxData=qammod(reshape(TxDataBits,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');%data=1*2560
Tx = Generate_2D_data_grid(N,M,TxData,data_grid);
Tx_Symb=Tx_addPilotSync(Tx,PilotBits,SyncBits,N,M_mod); 

%% OTSM modulation%%%%
Tx_tilda=Tx_Symb*Wn;              %equation (6) in [R1]   %Tx=X
%TxDataOtsmSymb=OtsmSignalModulation(Tx_tilda, NumFFT, NumCP);
tx_Data_signal=reshape(Tx_tilda,[],1);  %equation (7) in [R1]
tx_signal = [ ...
    tx_Data_signal(N*M-NumCP+1:N*M,1)
    tx_Data_signal];
flt1=rcosine(1,upsample,'fir/sqrt',0.05,64);%pulse shaper 
tx_signal2=rcosflt(tx_signal,1,upsample, 'filter', flt1); %3040(TxSignal)*4(upsample)+(513(flt1)-1)：因為捲積所以要-1