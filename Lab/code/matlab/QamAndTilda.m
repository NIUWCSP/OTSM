function XSymb_tilda=QamAndTilda(X,M_mod,M_bits,N,M,Wn)

QamXBits=qammod(reshape(X,M_bits,size(X,2)/2), M_mod,'gray','InputType','bit');
%QamSyncBits=reshape(QamSyncBits,[],1);
QamXBits=reshape(QamXBits,sqrt(size(QamXBits,2)),[]);%切成方形矩陣
%%加入只有同步的網格並做WHT
XGrid=zeros(N,M);
XGrid(1:size(QamXBits,1),1:size(QamXBits,2))=QamXBits;
XGrid_tilda=XGrid*Wn;
XSymb_tilda=reshape(XGrid_tilda(1:size(QamXBits,1),1:size(QamXBits,2)),[],1);

end