function tx_signal2 = Transmitter(upsample)
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

% Generate data symbols
global TxDataBits;
TxDataBits = randi([0,1],N_syms_perfram*M_bits,1);%TX的data
TxDataOtsmSymbMtx = qammod(reshape(TxDataBits,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit'); 
Tx = Generate_2D_data_grid(N,M,TxDataOtsmSymbMtx,data_grid);
TxDataOtsmSymb = reshape(TxDataOtsmSymbMtx, [], 1);

%% OTSM modulation%%%%
Tx_tilda=Tx*Wn;              %equation (6) in [R1]   %Tx=X
tx_signal2=reshape(Tx_tilda,N*M,1);  %equation (7) in [R1]