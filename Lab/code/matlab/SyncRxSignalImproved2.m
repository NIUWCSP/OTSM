function startIdx = SyncRxSignalImproved2(rxFrame,M_mod,N,M)

%% Set start index to an invalid number
startIdx = -1;
% M_mod: size of QAM constellation
M_bits = log2(M_mod);

%% Normalized WHT matrix
Wn=fwht(eye(N));  % Generate the WHT matrix
Wn=Wn./norm(Wn);  % normalize the WHT matrix

%% Construct the syncSig to be used for fine tuning

%%調變同步資料
SyncBits = GetSyncBits();%確保正確解讀接收到的數據
QamSyncBits=qammod(reshape(SyncBits,M_bits,size(SyncBits,2)/2), M_mod,'gray','InputType','bit');
QamSyncBits=reshape(QamSyncBits,sqrt(size(QamSyncBits,2)),[]);%切成方形矩陣
%%加入只有同步的網格並做WHT
SyncGrid=zeros(N,M);
SyncGrid(1:size(QamSyncBits,1),1:size(QamSyncBits,2))=QamSyncBits;
SyncGrid_tilda=SyncGrid*Wn;
SyncSymb_tilda=reshape(SyncGrid_tilda(1:size(QamSyncBits,1),1:size(QamSyncBits,2)),[],1);

syncSig = SyncSymb_tilda;
% 计算自相关
autoCorr = abs(xcorr(rxFrame, syncSig)); % 计算接收信号与前导信号的自相关

% 寻找峰值
[~, locs] = max(autoCorr); 

% 取最高峰值对应的位置作为帧起始位置
startIdx = locs + length(syncSig);