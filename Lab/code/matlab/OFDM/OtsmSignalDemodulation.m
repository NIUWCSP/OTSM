
function OfdmSymb = OtsmSignalDemodulation(ModOtsmSymbWithCP, ...
                                NumFFT, NumCP, NumDataCarrier,M_mod)

ModOtsmSymb = ModOtsmSymbWithCP(NumCP+1:end, :);

OfdmSymbFFT = fft(ModOtsmSymb) / sqrt(NumFFT);
%%在FFT的過程中，訊號的振幅會被FFT的點數所縮放。 為了在頻域中正確表示原始訊號的幅度

OtsmSymbFFT = qamdemod(OfdmSymbFFT,M_mod,'gray','OutputType','bit');

DemodOtsmSymb = [ ...
    OtsmSymbFFT(NumFFT/2+1:end, :);
    OtsmSymbFFT(1:NumFFT/2, :)];
OfdmSymb = [ ...
    DemodOtsmSymb(NumFFT/2-NumDataCarrier/2+1:NumFFT/2, :);
    DemodOtsmSymb(NumFFT/2+2:NumFFT/2+NumDataCarrier/2+1, :)];

