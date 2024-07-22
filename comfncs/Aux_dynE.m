
% ----------------------------------------------------------------------- %
% two-dimensional entry vehicle dynamic equations
% based on planet-centered rotating frame and north-east-down frame
% e = mu/r - 0.5V^2 is the independent variable
% ----------------------------------------------------------------------- %

function Xdot = Aux_dynE(e, X, sigs, e0, auxdata)

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
BAS = auxdata.BAS;

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
if BAS == 1 % linear
    sigf  = auxdata.sigf;
    sig0  = sigs(1);
    sigma = sig0 + (sigf - sig0)/(ef - e0)*(e - e0);

elseif BAS == 2 % exponential
    KEF   = auxdata.KEF;
    sigf  = auxdata.sigf;
    sig0  = sigs(1);
    sigma = exp(-(e-e0)/(ef-e0)*KEF)*(sig0-sigf) + sigf;

elseif BAS == 3 % exponential
    KEF   = auxdata.KEF;
    sig0  = sigs(1);
    sigma = sig0*exp(-(e-e0)/(ef-e0)*KEF);

elseif BAS == 4 % logistic
    KLF   = auxdata.KLF;
    sig0  = sigs(1);
    sigma = 2*sig0./(1 + exp(((e-e0)/(ef-e0)*KLF)));

elseif BAS == 5 % two-segment linear
    sigf = auxdata.sigf;
    sig0 = abs(sigs(1));
    sigm = abs(sigs(2));
    emid = 0.5*(e0+ef);
    if e < emid
        sigma = sig0 + (sig0 - sigm)/(e0 - emid)*(e - e0);
    else
        sigma = sigm + (sigf - sigm)/(ef - emid)*(e - emid);
    end

elseif BAS == 6 % quadratic
    sigf = auxdata.sigf;
    sig0 = abs(sigs(1));
    sigm = abs(sigs(2));

    a = 2*(sig0 - 2*sigm + sigf)/(e0-ef)^2;
    b = -(sig0*e0 + 3*sig0*ef - 4*sigm*e0 - 4*sigm*ef + 3*sigf*e0 ...
        + sigf*ef)/(e0-ef)^2;
    c = (sig0*e0*ef + sig0*ef^2 - 4*sigm*e0*ef + sigf*e0^2 ...
        + sigf*e0*ef)/(e0-ef)^2;
    sigma = e.^2.*a + b.*e + c;

elseif BAS == 7 % automatic linear
    sig0 = abs(sigs(1));
    sigf = abs(sigs(2));

    sigma = sig0 + (sigf - sig0)/(ef - e0)*(e - e0);
end

% Energy-based longitudinal entry dynamics
dedt = D*V;
dsdt = V*cos(gamma)/r;
drdt = V*sin(gamma);
dgdt = (L/V)*cos(sigma) + (V/r - g/V)*cos(gamma);

dsde = dsdt/dedt;
drde = drdt/dedt;
dgde = dgdt/dedt;

Xdot = [ dsde drde dgde]';
