function X=ICompeteISI(OtsmSymbWithCP,NumDataCarrier,M,NumCP)

OtsmSymb = OtsmSymbWithCP(NumCP+1:end, :);

OtsmSymbFFT = fft(OtsmSymb) / sqrt(M*2);
Xsize=length(OtsmSymbFFT);
%%在FFT的過程中，訊號的振幅會被FFT的點數所縮放。 為了在頻域中正確表示原始訊號的幅度
DemodOtsmSymb = [ ...
     OtsmSymbFFT(Xsize/2+1:end, :);
     OtsmSymbFFT(1:Xsize/2, :)];
%%將FFT輸出的頻譜從零頻率移到中心
X = [ ...
    DemodOtsmSymb(Xsize/2-NumDataCarrier/2+1:Xsize/2,:);
    DemodOtsmSymb(Xsize/2+1:Xsize/2+NumDataCarrier/2,:)];
%%選擇感興趣的資料子載波部分