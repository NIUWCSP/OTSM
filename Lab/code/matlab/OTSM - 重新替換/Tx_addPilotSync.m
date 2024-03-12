function TxRadioFrame=Tx_addPilotSync(Tx,PilotBits,SyncBits,N,M_mod)
%陳昱升
% M_mod: size of QAM constellation
M_bits = log2(M_mod);
NumDataN=N;

%%算出傳送端資料部分X相位

    for idx=1:N
        if Tx(idx,:)==zeros(1,N)
            NumDataN=NumDataN-1;
            %fprintf('%d:%d\n',idx,NumDataN);
        end
    end
%%Pilot QAMmodulation
QamPilotBits=qammod(reshape(PilotBits,M_bits,size(PilotBits,2)/2), M_mod,'gray','InputType','bit');
QamPilotBits=reshape(QamPilotBits,sqrt(size(QamPilotBits,2)),[]);%切成方形矩陣

%%Sync QAMmodulation
%QamSyncBits=qammod(reshape(SyncBits,M_bits,size(SyncBits,2)/2), M_mod,'gray','InputType','bit');
%QamSyncBits=reshape(QamSyncBits,sqrt(size(QamSyncBits,2)),[]);%切成方形矩陣

%%再資料後面加上 空白+Pilot+空白
TxRadioFrame=Tx;
%TxRadioFrame(NumDataN+1:NumDataN+size(QamSyncBits,2),1:size(QamSyncBits,2))=QamSyncBits;
TxRadioFrame(NumDataN+size(QamPilotBits,2)+1:NumDataN+size(QamPilotBits,2)*2,1:size(QamPilotBits,2))=QamPilotBits;

%for idy=1:N/size(SPS,2)
%TxRadioFrame(NumDataN+1:N,(idy-1)*size(SPS,2)+1:idy*size(SPS,2))=SPS;%在底下空白處加入SPS
%end
