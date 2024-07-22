% Mars_human_A.m
% ----------------------------------------------------------------------- %
% problem formulation for the atmospheric entry phase of the
% human Mars EDL mission (Mid lift-to-drag ratio rigid vehicle (MRV))
% ----------------------------------------------------------------------- %

% constants
rp = 3397e3;    % Mars radius, m
mu = 4.2828e13; % gravitational parameter, m^3/s^2
w  =  7.0882e-5; % spin rate, rad/s
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
h0 = 125000;  % m
r0  = rp + h0; % radius, m
theta0 =  -1.63*d2r;  % longitude, deg
phi0   =  -22.18*d2r; % latitude, deg
V0     = 4700;       % velocity, m/s
gamma0 =   -10.2*d2r;  % flight path angle in deg
psi0   =  4.6*d2r;    % heading angle in deg

% coordinate tranformation, inertial frame to rotating frame
[r0, theta0, phi0, V0, gamma0, psi0] = ...
    Aux_pci2pcr(r0, theta0, phi0, V0, gamma0, psi0, w);

% final condition (AIAA 2017-1898)
hf       = 2480; %3200;
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