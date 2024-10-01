
% ----------------------------------------------------------------------- %
% calculate heading angle error, which is an included angle between
% velocity vector projected onto local horizontal plane and
% target pointing vector projected onto local horizontal plane
% at current location
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %
% input  = state variables of vehicle and target location in rad
% output = heading angle error in rad

function del_psi = cal_delpsi(X, thetat, phit)

% state variable
r = X(1); theta = X(2); phi = X(3); V = X(4); gamma = X(5); psi = X(6);

% sin/cos of angles
cgam = cos(gamma); sgam = sin(gamma); cpsi = cos(psi); spsi = sin(psi);
cth  = cos(theta); sth  = sin(theta); cphi = cos(phi); sphi = sin(phi);
ctht = cos(thetat); stht = sin(thetat); cphit = cos(phit); sphit = sin(phit);

% ------ projection of velocity vector onto local horizontal plane ------ %
% 1. get velocity vector with given MCR component (NED frame)
% 2. projection onto the local-horizontal frame (= North East plane)
%    then third component = 0
% 3. express the velocity vector in planet-fixed frame (ex. ECEF)
% 4. normalize the velocity vector
Vp = V*cgam*[ -cth*sphi*cpsi - sth*spsi ; -sth*sphi*cpsi + cth*spsi ; cphi*cpsi ];
UVp = Vp/norm(Vp);

% position vector
r_nav = [ cth*cphi; sth*cphi; sphi];
Ur = r_nav/norm(r_nav);

% target vector on ground
r_tgt = [ ctht*cphit ; stht*cphit ; sphit ];
Del = r_tgt - r_nav;

% project Del onto horizontal plane
Delp = cross(cross(Ur,Del),Ur);
UDelp = Delp/norm(Delp);

% heading error magnitude
del_psi = acos(dot(UVp,UDelp));

% heading error direction
% if Vp is on the right side of Delp
% --> cross(UVp,UDelp) // Ur   --> + dot product
% if Vp is on the left side of Delp
% --> cross(UVp,UDelp) // -Ur   --> - dot product
sign_ck = sign(dot( cross(UVp,UDelp), Ur ));

if sign_ck < 0
    del_psi = -abs(del_psi);
end

del_psi = real(del_psi);