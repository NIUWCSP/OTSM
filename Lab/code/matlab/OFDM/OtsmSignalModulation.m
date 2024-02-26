function OtsmSymbWithCP = OtsmSignalModulation(Bits, NumFFT, NumCP, M_mod)
% Symbol mapping: QAM

NumOtsmSymb = size(Bits, 2);
MapSymb = zeros(size(qammod(Bits, M_mod,'gray','InputType','bit')));
for idx = 1:NumOtsmSymb
    MapSymb(:, idx) = qammod(Bits(:,idx), M_mod,'gray','InputType','bit');
end

%對抗ISI
MapFft = [ ...
    zeros(NumFFT/2-size(MapSymb,1)/2, NumOtsmSymb);%16個0
    MapSymb(1:end/2, :);
    zeros(1, NumOtsmSymb);%1個0
    MapSymb(end/2+1:end, :);
    zeros(NumFFT/2-size(MapSymb,1)/2-1, NumOtsmSymb)];%15個0

%前後一半對調
ReMapFft = [ ...
    MapFft(NumFFT/2+1:end, :);
    MapFft(1:NumFFT/2, :)];

%IFFT
OtsmSymb = ifft(ReMapFft) * sqrt(NumFFT); %"* sqrt(NumFFT)"取normalization

OtsmSymbWithCP = [ ...
    OtsmSymb(NumFFT-NumCP+1:end, :);
    OtsmSymb];
