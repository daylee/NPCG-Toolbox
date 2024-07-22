
% ----------------------------------------------------------------------- %
% calculate downrange and crossrange of entry vehicle to a target
% using the parameterized bank angle of the NPCG algorithm
% ----------------------------------------------------------------------- %
% input  = current state
% output = downrange and crossrange

function [DR, CR] = cal_z_3d(X0, e0, convar, BR, auxdata)

% auxdata
theta0   = auxdata.theta0;
phi0     = auxdata.phi0;
thetatgt = auxdata.thetatgt;
phitgt   = auxdata.phitgt;
ef       = auxdata.ef;
opt      = auxdata.opt;

% propagation
espan = [e0  ef];
[~,X] = ode45(@(e,X) Aux_dynE_3d(e, X, convar, BR, e0, auxdata),...
    espan, X0, opt);

% downrange and crossrange
theta = X(end,2);
phi = X(end,3);
[DR, CR] = cal_drdc(theta0, phi0, thetatgt, phitgt, theta, phi);