function OtsmSymbWithCP = OtsmSignalModulation(Bits, NumFFT, NumCP)
% Symbol mapping: QAM

NumOtsmSymb = size(Bits, 2);
Bits_Signal=reshape(Bits,[],1);

%對抗ISI
MapFft = [ ...
    zeros(NumFFT/2-size(Bits,1)/2, NumOtsmSymb);%16個0
    Bits(1:end/2, :);
    zeros(1, NumOtsmSymb);%1個0
    Bits(end/2+1:end, :);
    zeros(NumFFT/2-size(Bits,1)/2-1, NumOtsmSymb)];%15個0

%前後一半對調
ReMapFft = [ ...
    MapFft(NumFFT/2+1:end, :);
    MapFft(1:NumFFT/2, :)];

%IFFT
OtsmSymb = ifft(ReMapFft) * sqrt(NumFFT); %"* sqrt(NumFFT)"取normalization

OtsmSymbWithCP = [ ...
    OtsmSymb(NumFFT-NumCP+1:end, :);
    OtsmSymb];
