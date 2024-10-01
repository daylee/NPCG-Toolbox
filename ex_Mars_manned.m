
% ----------------------------------------------------------------------- %
% problem formulation for the atmospheric entry phase of the
% human Mars EDL mission (Mid lift-to-drag ratio rigid vehicle (MRV))
% Main reference: Jiang et al., Acta Astronautica 163 (2019) 114-129
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %

% constants
rp = 3397e3;    % Mars radius, m
mu = 4.2828e13; % gravitational parameter, m^3/s^2
w  = 7.0882e-5; % spin rate, rad/s
d2r = pi/180;
r2d = 180/pi;

% planet / 1: Earth, 2: Mars
pn = 2;

% vehicle parameters, L/D = 0.54, B = 379 kg/m^2
S = 160;     % cross-sectional area, m^2
m = 58800;   % mass, kg
CL = 0.5236; % lift cofficient
CD = 0.9696; % drag cofficient

% entry interface in planet-centered-inertial frame (Jiang et al., 2019)
h0     = 125000;  % m
r0     = rp + h0; % radius, m
theta0 =  -1.63*d2r;  % longitude, deg
phi0   =  -22.18*d2r; % latitude, deg
V0     = 4700;       % velocity, m/s
gamma0 =   -10.2*d2r;  % flight path angle in deg
psi0   =  4.6*d2r;    % heading angle in deg

% coordinate tranformation, inertial frame to rotating frame
[r0, theta0, phi0, V0, gamma0, psi0] = ...
    cal_pci2pcr(r0, theta0, phi0, V0, gamma0, psi0, w);

% final condition (AIAA 2017-1898)
hf       = 2480;
rf       = hf + rp;
Vf       =  480;
thetatgt =   0*d2r;
phitgt   =  0*d2r;

% path constraints (AAS 19-619)
Amax    = 4;   % g-load, g
qmax    = 13;  % dynamic pressure, kPa
Qdotmax = 500; % heating rate, kW/m^2

% heating rate constant
kq = 5.3697e-8;
kqN = 0.5;
kqM = 3.15;

% bank angle constraints
BAL = 1;
siglmt = 180*d2r;   % bank magnitude limit
sigdlmt =  15*d2r;  % bank rate limit
sigddlmt = 15*d2r;  % bank acc limit

% guidance activation time
GAT= 170;

% bank reversal logic
BRL = 2;
dlpsT = [6 3]*d2r;
KBR = 7;

% bank angle parameterization options
if BAP == 1

    sig0 = 70*d2r; % initial guess of the parameter to be found
    sigf = 40*d2r; % design parameter

    sigs = sig0;
    auxdata.sigf = sigf;

elseif BAP == 2

    sig0 = 70*d2r; % initial guess of the parameter to be found
    sigf = 20*d2r; % design parameter
    KEF  = 1.2;    % design parameter

    sigs = sig0;
    auxdata.sigf = sigf;
    auxdata.KEF = KEF;

elseif BAP == 3

    sig0 = 70*d2r; % initial guess of the parameter to be found

    KEF  = 0.9; % design parameter

    sigs = sig0;
    auxdata.KEF = KEF;

elseif BAP == 4

    sig0 = 70*d2r; % initial guess of the parameter to be found

    KLF  = 1.3; % design parameter

    sigs = sig0;
    auxdata.KLF = KLF;
    
end