
% ----------------------------------------------------------------------- %
% calculate range-to-go at the current state
% using the parameterized bank angle of the NPCG algorithm
% ----------------------------------------------------------------------- %
% input  = current state
% output = range-to-go in rad

function z = cal_G(X0, convar, e0, auxdata)

% auxdata
thetatgt = auxdata.thetatgt;
phitgt   = auxdata.phitgt;
ef       = auxdata.ef;
opt      = auxdata.opt;

% state vector
r0     = X0(1);
theta0 = X0(2);
phi0   = X0(3);
gamma0 = X0(5);

% current range-to-go
sf = cal_sphdis(theta0, phi0, thetatgt, phitgt);

% initial condition
x0 = [0 r0 gamma0];

% energy span
espan = [e0  ef];

% propagation
[~,X] = ode45(@(e,X) Aux_dynE(e, X, convar, e0, auxdata), espan, x0, opt);

sef = X(end,1) - sf;

z = sef;