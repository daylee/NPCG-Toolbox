
% ----------------------------------------------------------------------- %
% problem formulation for the atmospheric entry phase of the
% human Earth EDL mission, Apollo 10 command module
% Main reference: Szelc 1969, NASA TM-X-69555
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %

% constants
rp = 6378.137e3;      % Earth radius, m
mu = 398600.435507e9; % gravitational parameter, m^3/s^2
w  = 2*pi/86400;      % spin rate, rad/s
d2r = pi/180;
r2d = 180/pi;

% planet / 1: Earth, 2: Mars
pn = 1;

% vehicle parameters
S = 12.017;   % cross-sectional area, m^2
m = 5498.22;  % mass, kg
CL = 0.40815; % lift cofficient
CD = 1.2569;  % drag cofficient

% entry interface in planet-centered-inertial frame (Szelc 1969)
r0     = 6498.270e3;     % geocentric radius, m
theta0 = 174.24384*d2r;  % longitude, deg
phi0   = -23.51457*d2r;  % geocentric latitude, m
V0     = 11.06715e3;     % inertial velocity, m/s
gamma0 = -6.6198381*d2r; % inertial flight path angle
psi0   = 71.9317*d2r;    % inertial heading angle

% coordinate tranformation, inertial frame to rotating frame
[r0, theta0, phi0, V0, gamma0, psi0] = ...
    cal_pci2pcr(r0, theta0, phi0, V0, gamma0, psi0, w);

% final condition
hf       = 5683;
rf       = hf + rp;
Vf       = 119.67;
thetatgt = 194.8501*d2r;
phitgt   = -15.4121*d2r;

% path constraints (Szelc 1969, Pavlosky and Leger 1974)
Amax    = 10;   % g-load, g
qmax    = 30;   % dynamic pressure, kPa
Qdotmax = 6000; % heating rate, kW/m^2

% heating rate constant (Young and Smith 1967)
BTU2kW = 11.356539; % 1 BTU/s.ft2 to kW/m^2
kq = BTU2kW*20*1e-9;
kqN = 0.5;
kqM = 3;

% bank angle constraints
BAL = 1;
siglmt = 180*d2r;   % bank magnitude limit
sigdlmt = 15*d2r;  % bank rate limit
sigddlmt = 15*d2r;  % bank acc limit

% guidance activation time
GAT= 90;

% bank reversal logic
BRL = 1;
dlpsT = [5 7]*d2r;
KBR = 6;

% bank angle parameterization options
if BAP == 1

    sig0 = 70*d2r; % initial guess of the parameter to be found
    sigf = 30*d2r; % design parameter

    sigs = sig0;
    auxdata.sigf = sigf;

elseif BAP == 2

    sig0 = 70*d2r; % initial guess of the parameter to be found
    sigf = 10*d2r; % design parameter
    KEF  = 0.9;    % design parameter

    sigs = sig0;
    auxdata.sigf = sigf;
    auxdata.KEF = KEF;

elseif BAP == 3

    sig0 = 70*d2r; % initial guess of the parameter to be found
    KEF  = 1.3; % design parameter

    sigs = sig0;
    auxdata.KEF = KEF;

elseif BAP == 4

    sig0 = 70*d2r; % initial guess of the parameter to be found

    KLF  = 1.3; % design parameter

    sigs = sig0;
    auxdata.KLF = KLF;

end