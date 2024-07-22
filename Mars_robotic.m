
% ----------------------------------------------------------------------- %
% problem formulation for the atmospheric entry phase of the
% robotic Mars EDL mission (Mars science laboratory)
% ----------------------------------------------------------------------- %

% constants
rp = 3397e3;    % Mars radius, m
mu = 4.2828e13; % gravitational parameter, m^3/s^2
w  = 7.0882e-5; % spin rate, rad/s
d2r = pi/180;
r2d = 180/pi;

% planet / 1: Earth, 2: Mars
pn = 2;

% vehicle parameters, L/D = 0.24, B = 146 kg/m^2
S = 15.588; % cross-sectional area, m^2
m = 3300;   % mass, kg
CL = 0.348; % lift cofficient
CD = 1.45;  % drag cofficient

% entry interface in planet-centered-inertial frame (Dutta and Braun 2014)
r0     = rp + 125e3;   % radius, m
theta0 = 126.72*d2r;   % longitude, deg
phi0   = -3.9186*d2r;  % latitude, deg
V0     = 6083.3;       % velocity, m/s
gamma0 = -15.4892*d2r; % flight path angle in deg
psi0   = 93.2065*d2r;  % heading angle in deg

% coordinate tranformation, inertial frame to rotating frame
[r0, theta0, phi0, V0, gamma0, psi0] = ...
    Aux_pci2pcr(r0, theta0, phi0, V0, gamma0, psi0, w);

% final condition 
hf       = 12000;
rf       = hf + rp;
Vf       = 406;
thetatgt = 137.431*d2r;
phitgt   = -4.552*d2r;

% path constraints limits
Amax    = 15;   % g-load, g
qmax    = 21;   % dynamic pressure, kPa
Qdotmax = 1000; % heating rate, kW/m^2

% heating rate constant
kq = 5.3697e-8;
kqN = 0.5;
kqM = 3.15;