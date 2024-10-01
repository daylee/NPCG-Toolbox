
% ----------------------------------------------------------------------- %
% 3-degree of freedom equations of motion of the entry vehicle
% based on planet-centered rotating frame and north-east-down frame
% time is the independent variable
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %

function Xdot = aux_dyn(t, X, sigma, auxdata)

% auxdata
pn = auxdata.pn;
rp = auxdata.rp;
mu = auxdata.mu;
w  = auxdata.w;
S  = auxdata.S;
m  = auxdata.m;
CL = auxdata.CL;
CD = auxdata.CD;
DU = auxdata.DU;
VU = auxdata.VU;
TU = auxdata.TU;

% state
r = X(1); theta = X(2); phi = X(3);
V = X(4); gamma = X(5); psi = X(6);

% aerodynamic forces and gravity
Vreal = V*VU; % m/s
hreal = (r - rp)*DU; % m
rho = cal_airdens(hreal, pn); % kg/m^3
q = 0.5*rho*Vreal^2;
D = q*S*CD/m;
L = q*S*CL/m;
D = D/DU*TU^2; % unitless drag
L = L/DU*TU^2; % unitless lift
g = mu*r^-2;   % unitless gravity

% time domain 3-DoF entry dynamics
rdot     = V*sin(gamma);
thetadot = V*cos(gamma)*sin(psi)/(r*cos(phi));
phidot   = (V/r)*cos(gamma)*cos(psi);
Vdot   = - D - g*sin(gamma)...
    + w^2*r*cos(phi)*(sin(gamma)*cos(phi) - cos(gamma)*sin(phi)*cos(psi));
gamdot = (L/V)*cos(sigma) + (V/r - g/V)*cos(gamma)...
    + 2*w*cos(phi)*sin(psi) ...
    + w^2*r/V*cos(phi)*(cos(gamma)*cos(phi) + sin(gamma)*sin(phi)*cos(psi));
psidot = (L*sin(sigma))/(V*cos(gamma)) + (V/r)*sin(psi)*cos(gamma)*tan(phi)...
    - 2*w*(tan(gamma)*cos(phi)*cos(psi) - sin(phi)) ...
    + w^2*r/(V*cos(gamma))*sin(phi)*cos(phi)*sin(psi);

Xdot = [ rdot  thetadot  phidot  Vdot  gamdot  psidot ]';