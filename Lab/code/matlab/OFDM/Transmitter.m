function tx_signal2 = Transmitter(upsample)
NumFFT = 128;%V3 FFT轉換的點數
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
% Generate synchronization symbols
SyncBits = GetSyncBits();%Preamble的data
SyncOfdmSymb = OfdmSignalModulation(SyncBits, NumFFT, 0);

% Generate pilot symbols
PilotBits = GetPilotBits();%Preamble的data
TxPilotOfdmSymb = OfdmSignalModulation(PilotBits, NumFFT, 0);

% Generate data symbols
global TxDataBits;
TxDataBits = randi([0,1],N_syms_perfram*M_bits,1);%TX的data
TxDataOtsmSymbMtx = qammod(reshape(TxDataBits,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit'); 
TxDataOtsmSymb = reshape(TxDataOtsmSymbMtx, [], 1);

% Reconstruct transmission signal
TxSignal = [ ...
    SyncOfdmSymb(1:NumSyncPreamble);%"SyncOfdmSymb"三次目的是要方便同步(去除phase offset)
    SyncOfdmSymb(1:NumSyncPreamble);
    SyncOfdmSymb;
    TxPilotOfdmSymb;%通道估測
    TxPilotOfdmSymb;
    TxDataOtsmSymb];
TxDataBits = Generate_2D_data_grid(N,M,TxDataOtsmSymbMtx,data_grid);
%% OTSM demodulation%%%%
tx_signal2=TxDataBits*Wn