function BER ( TxDataBits , RxDataBits , NoFoundDataTimes , AllFoundDataTimes)
    
    TxDataBitsReshape = reshape(TxDataBits, [], 1);
    RxDataBitsReshape = reshape(RxDataBits, [], 1);
    BER = sum(abs(TxDataBitsReshape-RxDataBitsReshape));
    ER = BER/length(TxDataBitsReshape);
    subplot(235);
    axis off;
    
    text(1.5,1.0,['BER: ',num2str(ER*100),' %']);
    text(1.5,0.8,['ER: ',num2str(BER),' bits']);
    text(1.5,0.6,['NoFoundData: ',num2str(NoFoundDataTimes),' 次']);
    text(1.5,0.4,['AllFoundData: ',num2str(AllFoundDataTimes),' 次']);
    text(1.5,0.2,['無收訊率: ',num2str(NoFoundDataTimes / AllFoundDataTimes *100),' %']);
    global BERData;
    BERData = ER;
end
