
% ----------------------------------------------------------------------- %
% predictive lateral guidance based on
% Smith, K. M. (2016, February).
% Predictive lateral logic for numerical entry guidance algorithms.
% In AAS/AIAA Space Flight Mechanics Meeting (No. JSC-CN-35110-1).
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %
% input  = current state, previous bank sign,
%          parameterzied bank angle of the NPCG algorithm
% output = bank sign

function BR = cal_BR_prdt(X0, e0, bank, BRprev, auxdata)

% auxdata
KBR = auxdata.KBR;

% crossrange error based on the current bank profile (mag + sign)
[~, CR1] = cal_z_3d(X0, e0, bank, BRprev, auxdata);

% crossrange error based on the opposite sign of bank profile
[~, CR2] = cal_z_3d(X0, e0, bank, BRprev*-1, auxdata);

ratio = abs(CR1/CR2);

if ratio < KBR
    BR = BRprev;
else
    BR = BRprev*-1;
end