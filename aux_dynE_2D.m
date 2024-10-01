
% ----------------------------------------------------------------------- %
% 2-degree of freedom equations of motion of the entry vehicle
% based on planet-centered rotating frame and north-east-down frame
% energy-like variable, e = mu/r - 0.5V^2, is the independent variable
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %

function Xdot = aux_dynE_2D(e, X, sigs, e0, auxdata)

% auxdata
pn  = auxdata.pn;
rp  = auxdata.rp;
mu  = auxdata.mu;
S   = auxdata.S;
m   = auxdata.m;
CL  = auxdata.CL;
CD  = auxdata.CD;
DU  = auxdata.DU;
VU  = auxdata.VU;
TU  = auxdata.TU;
ef  = auxdata.ef;
BAP = auxdata.BAP;

% state
s = X(1); r = X(2); gamma = X(3); V = sqrt(2*(mu/r - e));

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

% bank angle parameterization
if BAP == 1 % linear
    sigf  = auxdata.sigf;
    sig0  = sigs(1);
    sigma = sig0 + (sigf - sig0)/(ef - e0)*(e - e0);

elseif BAP == 2 % exponential 1
    KEF   = auxdata.KEF;
    sigf  = auxdata.sigf;
    sig0  = sigs(1);
    sigma = exp(-(e-e0)/(ef-e0)*KEF)*(sig0-sigf) + sigf;

elseif BAP == 3 % exponential 2
    KEF   = auxdata.KEF;
    sig0  = sigs(1);
    sigma = sig0*exp(-(e-e0)/(ef-e0)*KEF);

elseif BAP == 4 % logistic
    KLF   = auxdata.KLF;
    sig0  = sigs(1);
    sigma = 2*sig0./(1 + exp(((e-e0)/(ef-e0)*KLF)));

end

% energy domain 2-DoF entry dynamics
dedt = D*V;
dsdt = V*cos(gamma)/r;
drdt = V*sin(gamma);
dgdt = (L/V)*cos(sigma) + (V/r - g/V)*cos(gamma);

dsde = dsdt/dedt;
drde = drdt/dedt;
dgde = dgdt/dedt;

Xdot = [ dsde drde dgde]';