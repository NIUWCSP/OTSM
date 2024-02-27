function tx_signal2 = Transmitter(upsample)
NumFFT = 64;%V3 FFT轉換的點數
NumSyncPreamble = 32;%V3 同步的前綴，Preamble：防干擾+同步+通道估測(已知的頻域資料)
NumCP = 16;%V3 CP：循環前綴(NumFFT = 128後16貼回前面)，CP：避免ISI(多路徑干擾)(未知的時域訊號)

%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);

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

%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix

%% Transmitter
% Generate pilot symbols
PilotBits = GetPilotBits();%Preamble的data
TxPilotOtsmSymb = OtsmSignalModulation(PilotBits, NumFFT, 0 , M_mod);

% Generate synchronization symbols
SyncBits = GetSyncBits();%Preamble的data
SyncOtsmSymb = OtsmSignalModulation(SyncBits, NumFFT, 0 , M_mod);

% Generate data symbols
global TxDataBits;
TxDataBits = randi([0,1],N_syms_perfram,1);%TX的data
TxDataBitsGrid=reshape(TxDataBits,N,[]);
TxDataOtsm = OtsmSignalModulation(TxDataBitsGrid, NumFFT, 0 , M_mod);
Tx = Generate_2D_data_grid(N,M,TxDataOtsm,data_grid);
TxDataOtsmSymb = reshape(Tx, [], 1);

%% OTSM modulation%%%%
Tx_tilda=Tx*Wn;              %equation (6) in [R1]   %Tx=X
Tx_tilda_Pilot=Tx_addPilot(Tx_tilda,TxPilotOtsmSymb);
tx_signal=reshape(Tx_tilda_Pilot,N*M,1);  %equation (7) in [R1]
tx_signal = [ ...
    SyncOtsmSymb(1:NumSyncPreamble);%"SyncOtsmSymb"三次目的是要方便同步
    SyncOtsmSymb(1:NumSyncPreamble);
    SyncOtsmSymb;
    TxPilotOtsmSymb;%通道估測
    TxPilotOtsmSymb;
    TxDataOtsmSymb];
flt1=rcosine(1,upsample,'fir/sqrt',0.05,64);%pulse shaper %1*513
tx_signal2=rcosflt(tx_signal,1,upsample, 'filter', flt1); %3040(TxSignal)*4(upsample)+(513(flt1)-1)：因為捲積所以要-1
