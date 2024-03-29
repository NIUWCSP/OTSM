function ReMapFft=QamAndTilda(X,M_mod,M_bits,N,M,Wn)

QamXBits=reshape(qammod(reshape(X,M_bits,size(X,2)/2), M_mod,'gray','InputType','bit'),[],1);
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
    QamXBits(1:end/2, :);
    zeros(1,1);%1個0
    QamXBits(end/2+1:end, :);
    zeros(9,1)];%9個0

%前後一半對調
ReMapFft = [ ...
    MapFft(size(MapFft,1)/2+1:end, :);
    MapFft(1:size(MapFft,1)/2, :)];
end