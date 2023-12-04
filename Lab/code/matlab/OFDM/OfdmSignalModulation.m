
function OfdmSymbWithCP = OfdmSignalModulation(Bits, NumFFT, NumCP)
% Symbol mapping: BPSK
BpskModObj = comm.BPSKModulator('PhaseOffset', pi/4);%'PhaseOffset', pi/4：0度偏移成45度、180度偏移成225度

NumOfdmSymb = size(Bits, 2);
MapSymb = zeros(size(Bits));
for idx = 1:NumOfdmSymb
    MapSymb(:, idx) = step(BpskModObj,Bits(:, idx));
end

%對抗ISI
MapFft = [ ...
    zeros(NumFFT/2-length(MapSymb)/2, NumOfdmSymb);%10個0
    MapSymb(1:end/2, :);
    zeros(1, NumOfdmSymb);%1個0
    MapSymb(end/2+1:end, :);
    zeros(NumFFT/2-length(MapSymb)/2-1, NumOfdmSymb)];%9個0

%前後一半對調
ReMapFft = [ ...
    MapFft(NumFFT/2+1:end, :);
    MapFft(1:NumFFT/2, :)];

%IFFT
OfdmSymb = ifft(ReMapFft) * sqrt(NumFFT); %"* sqrt(NumFFT)"取normalization

OfdmSymbWithCP = [ ...
    OfdmSymb(NumFFT-NumCP+1:end, :);
    OfdmSymb];
