function tx_signal2 = Transmitter(upsample)
NumFFT = 128;%V3 FFT轉換的點數
NumSyncPreamble = 32;%V3 同步的前綴，Preamble：防干擾+同步+通道估測(已知的頻域資料)
NumCP = 16;%V3 CP：循環前綴(NumFFT = 128後16貼回前面)，CP：避免ISI(多路徑干擾)(未知的時域訊號)

%% Transmitter
% Generate synchronization symbols
SyncBits = GetSyncBits();%Preamble的data
SyncOfdmSymb = OfdmSignalModulation(SyncBits, NumFFT, 0);

% Generate pilot symbols
PilotBits = GetPilotBits();%Preamble的data
TxPilotOfdmSymb = OfdmSignalModulation(PilotBits, NumFFT, 0);

% Generate data symbols
global TxDataBits;
TxDataBits = GetDataBits();%TX的data
TxDataOfdmSymbMtx = OfdmSignalModulation(TxDataBits, NumFFT, NumCP);
TxDataOfdmSymb = reshape(TxDataOfdmSymbMtx, [], 1);

% Reconstruct transmission signal
TxSignal = [ ...
    SyncOfdmSymb(1:NumSyncPreamble);%"SyncOfdmSymb"三次目的是要方便同步(去除phase offset)
    SyncOfdmSymb(1:NumSyncPreamble);
    SyncOfdmSymb;
    TxPilotOfdmSymb;%通道估測
    TxPilotOfdmSymb;
    TxDataOfdmSymb];
flt1=rcosine(1,upsample,'fir/sqrt',0.05,64);%pulse shaper %1*513
tx_signal2=rcosflt(TxSignal,1,upsample, 'filter', flt1); %3040(TxSignal)*4(upsample)+(513(flt1)-1)：因為捲積所以要-1