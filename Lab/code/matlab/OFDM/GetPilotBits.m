
function PilotBits = GetPilotBits(InvBits)

if nargin == 0
    InvBits = false;
end

% 陳
PilotBits = [...
    0 1 0 1 0 1 0 1 ...
    1 0 1 0 1 0 1 0 ...
    0 1 0 1 0 1 0 1 ...
    1 0 1 0 1 0 1 0 ...
    0 1 0 1 0 1 0 1 ...
    1 0 1 0 1 0 1 0 ...
    0 1 0 1 0 1 0 1 ...
    1 0 1 0 1 0 1 0 ...
    0 1 0 1 0 1 0 1 ...
    1 0 1 0 1 0 1 0 ...
    0 1 0 1 0 1 0 1 ...
    1 0 1 0 1 0 1 0 ...
    0 1 0 1 0 1 0 1 ...
    1 0 1 0 1 0 1 0 ...
    0 1 0 1 0 1 0 1 ...
    1 0 1 0 1 0 1 0
    ];
if InvBits
    PilotBits = ~PilotBits;
end