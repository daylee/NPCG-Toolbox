
% ----------------------------------------------------------------------- %
% coordinate transformation 
% from planet-centered-inertial (PCI) to planet-centered-rotating (PCR)
% based on north-east-down (NED) frame.
% heading angle (psi) is measured from North in clockwise
% ----------------------------------------------------------------------- %
% input  = state variables represented in PCI frame
% output = state variables represented in PCR frame

function [rn, thetan, phin, Vn, gamman, psin] = ...
    cal_pci2pcr(ri, thetai, phii, Vi, gammai,  psii, w)

% angles
sth  = sin(thetai); cth = cos(thetai);  sphi = sin(phii); cphi = cos(phii);
cgam = cos(gammai); sgam = sin(gammai); cpsi = cos(psii); spsi = sin(psii);

% Assumption!!! PCI and PCR are overlapped at the moment
thetan = thetai;
phin = phii;

% % radial vector in NED
% rN = [0;0;-ri];
% 
% % transformation from NED to PCI/PCR
% CPN = [-cth*sphi, -sth, -cphi*cth;
% -sphi*sth, cth, -cphi*sth;
% cphi, 0, -sphi];
% 
% % radial vector in PCI
% rI = CPN*rN;

% radial vector in PCI
rI = [ ri*cphi*cth;
    ri*cphi*sth;
    ri*sphi];

% % velocity vector in velocity-pointing frame
% vV = [Vi;0;0];
% 
% % transformation from velocity-pointing frame to PCI/PCR
% CPV = [cphi*cth*sgam-cgam*spsi*sth-cgam*cpsi*cth*sphi, cth*sphi*spsi-cpsi*sth, -sgam*spsi*sth-cgam*cphi*cth-cpsi*cth*sgam*sphi;
% cgam*cth*spsi+cphi*sgam*sth-cgam*cpsi*sphi*sth, cpsi*cth+sphi*spsi*sth, cth*sgam*spsi-cgam*cphi*sth-cpsi*sgam*sphi*sth;
% sgam*sphi+cgam*cphi*cpsi, -cphi*spsi, cphi*cpsi*sgam-cgam*sphi];
% 
% % velocity vector in PCI
% vI = CPV*vV;

% velocity vector in PCI
vI = [ Vi*( cth*cphi*sgam - sth*spsi*cgam - cth*sphi*cgam*cpsi );
    Vi*( sth*cphi*sgam + cth*spsi*cgam - sth*sphi*cgam*cpsi );
    Vi*( sphi*sgam + cphi*cgam*cpsi )];

% planet rotation vector in PCI
wI = w*[ 0 ; 0 ; 1 ];

% velocity vector in PCR (relative velocity)
vR = vI - cross(wI,rI);

% radius
rn = norm(rI);

% relative speed
Vn = norm(vR);

% planet-relative flight-path angle
gamman = asin( dot(rI,vR)/(rn*Vn) );

% planet-relative heading angle
sgamn = sin(gamman);
cgamn = cos(gamman);
psin = acos( (vR(3) - Vn*sgamn*sphi)/(Vn*cgamn*cphi));
