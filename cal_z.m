
% ----------------------------------------------------------------------- %
% calculate range-to-go using the parameterized bank angle
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %
% input  = current state
% output = range-to-go in rad

function z = cal_z(X0, sigs, e0, auxdata)

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
[~,X] = ode45(@(e,X) aux_dynE_2D(e, X, sigs, e0, auxdata), espan, x0, opt);

sef = X(end,1) - sf;

z = sef;