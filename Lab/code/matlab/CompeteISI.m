function OtsmSymbWithCP=CompeteISI(X,NumCP,N,M,Wn)

% QamXBits=reshape(qammod(reshape(X,M_bits,size(X,2)/2), M_mod,'gray','InputType','bit'),[],1);
%QamSyncBits=reshape(QamSyncBits,[],1);

%%加入只有同步的網格並做WHT
% QamXBitsGrid=reshape(QamXBits,sqrt(size(QamXBits,1)),[]);%切成方形矩陣
% XGrid=zeros(N,M);
% XGrid(1:size(QamXBitsGrid,1),1:size(QamXBitsGrid,2))=QamXBitsGrid;
% XGrid_tilda=XGrid*Wn;
% XSymb_tilda=reshape(XGrid_tilda(1:size(QamXBitsGrid,1),1:size(QamXBitsGrid,2)),[],1);

%對抗ISI
MapFft = [ ...
    zeros(10,1);%10個0
    X(1:end/2, :);
    zeros(1,1);%1個0
    X(end/2+1:end, :);
    zeros(9,1)];%9個0

%前後一半對調
ReMapFft = [ ...
    MapFft(size(MapFft,1)/2+1:end, :);
    MapFft(1:size(MapFft,1)/2, :)];
%%IFFT
%OtsmSymb =  ReMapFft * sqrt(M*2);
OtsmSymb = ifft(ReMapFft) * sqrt(M*3); %"* sqrt(M*2)"取normalization
OtsmSymbWithCP = [ ...
    OtsmSymb(size(OtsmSymb,1)-NumCP+1:end, :);
    OtsmSymb];

%%嘗試不用IFFT改用WHT
% ReMapFftGrid=reshape(ReMapFft,[],M);
% OtsmSymbGrid=zeros(N,M);
% OtsmSymbGrid(1:size(ReMapFftGrid,1),1:M)=ReMapFftGrid;
% OtsmSymbGrid_tilda=OtsmSymbGrid*Wn;
% FinalOtsmSymb=reshape(OtsmSymbGrid_tilda(1:size(ReMapFftGrid,1),1:M),[],1);
% OtsmSymbWithCP = [ ...
%     FinalOtsmSymb(size(FinalOtsmSymb,1)-NumCP+1:end, :);
%     FinalOtsmSymb];


end