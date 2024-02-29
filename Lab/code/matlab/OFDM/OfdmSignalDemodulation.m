
function OfdmSymb = OtsmSignalDemodulation(ModOtsmSymbWithCP, ...
                                NumFFT, NumCP, NumDataCarrier,M_mod,N_bits_perfram)

ModOtsmSymb = ModOtsmSymbWithCP(NumCP+1:end, :);

OfdmSymbFFT = reshape(qamdemod(ModOtsmSymb,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);
DemodOfdmSymb = [ ...
    OfdmSymbFFT(NumFFT/2+1:end, :);
    OfdmSymbFFT(1:NumFFT/2, :)];
OfdmSymb = [ ...
    DemodOfdmSymb(NumFFT/2-NumDataCarrier/2+1:NumFFT/2, :);
    DemodOfdmSymb(NumFFT/2+2:NumFFT/2+NumDataCarrier/2+1, :)];
