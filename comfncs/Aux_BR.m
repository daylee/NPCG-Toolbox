
% ----------------------------------------------------------------------- %
% bank reversal algorithm based on
% Tu, K. Y., Munir, M. S., Mease, K. D., & Bayard, D. S. (2000).
% Drag-based predictive tracking guidance for Mars precision landing.
% Journal of Guidance, Control, and Dynamics, 23(4), 620-628.
% ----------------------------------------------------------------------- %
% input  = current state, previous bank sign
% output = bank sign

function BR = Aux_BR(X, BR, ex_case, auxdata)

% auxdata
rp       = auxdata.rp;
DU       = auxdata.DU;
VU       = auxdata.VU;
thetatgt = auxdata.thetatgt;
phitgt   = auxdata.phitgt;

% state variables
r = X(1); theta = X(2); phi = X(3); V = X(4); gamma = X(5); psi = X(6);

% sin/cos of angles
cgam = cos(gamma); cpsi = cos(psi); cth  = cos(theta); cphi = cos(phi);
sgam = sin(gamma); spsi = sin(psi); sth  = sin(theta); sphi = sin(phi);

% ----------------------------------------------------------------------- %
%          velocity vector projection onto local horizontal plane
% ----------------------------------------------------------------------- %
% velocity vector projection onto NE plane --> NED to ECEF
Vp = V*cgam*[-cth*sphi*cpsi-sth*spsi; -sth*sphi*cpsi+cth*spsi; cphi*cpsi];
UVp = Vp/norm(Vp);

% position vector inertial frame
rnav = [ r*cos(theta)*cphi; r*sin(theta)*cphi; r*sphi];

% target vector inertial frame
rtgt = rp*[cos(thetatgt)*cos(phitgt);sin(thetatgt)*cos(phitgt);sin(phitgt)];

Del = rtgt - rnav;
Ur = rnav/norm(rnav);

% ----------------------------------------------------------------------- %
%       target pointing vector projection onto local horizontal plane
% ----------------------------------------------------------------------- %
Delp = cross(cross(Ur,Del),Ur);
UDelp = Delp/norm(Delp);

% heading angle error magnitude
delpsi = acos(dot(UVp,UDelp));

% if target is on the left side from the Vp
% --> cross(UVp,UDelp) // Ur --> + dot product
% if target is on the right side from the Vp
% --> cross(UVp,UDelp) // -Ur --> - dot product
sign_ck = dot( cross(UVp,UDelp), Ur );

% current range-to-go
R2go = cal_sphdis(theta, phi, thetatgt, phitgt)*DU*1e-3; % km

% current velocity
Vreal = V*VU;

d2r = pi/180;
% define heading angle error limit
if ex_case == 1 % Mars robotic

    if R2go > 100
        del_psi_lim = 4*d2r;
    else
        del_psi_lim = 1*d2r;
    end

elseif ex_case == 2 % Mars human

    if R2go > 150  % [km]
        del_psi_lim = 8*d2r;
    elseif R2go <= 150 && R2go > 50
        del_psi_lim = 3*d2r;
    elseif R2go <= 50
        del_psi_lim = 1*d2r;
    else
        return
    end

elseif ex_case == 3 % Apollo 10
    if R2go > 500
        del_psi_lim = 3*d2r;
    elseif R2go <= 500 && R2go > 10
        del_psi_lim = 14*d2r;
    elseif R2go <= 10
        del_psi_lim = 30*d2r;
    else
        return
    end

elseif ex_case == 4 % CAV-H mission

    if Vreal >= 6e3
        del_psi_lim = 10.5*d2r;
    elseif Vreal < 6e3 && Vreal >=5e3
        del_psi_lim = 18*d2r;
    elseif Vreal < 5e3 && Vreal >=3e3
        del_psi_lim = 20*d2r;
    else
        del_psi_lim = 15*d2r;
    end
end

% if current heading angle error exceeds the limit
if abs(delpsi) >= del_psi_lim
    if sign_ck > 0
        % target is on the left side of the flight direction
        % --> turn left --> negative roll
        BR = -1;
    else
        % target is on the right side of the flight direction
        % --> turn right --> positive roll
        % and negative heading angle erorr
        delpsi = -abs(delpsi);
        BR = 1;
    end
end
