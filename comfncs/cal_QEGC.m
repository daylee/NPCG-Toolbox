
% ----------------------------------------------------------------------- %
% calculate phugoid oscillation damping term based on
% Lu, P., Forbes, S., & Baldwin, M. (2013).
% Gliding guidance of high L/D hypersonic vehicles.
% In AIAA Guidance, Navigation, and Control (GNC) Conference (p. 4648).
% ----------------------------------------------------------------------- %
% input  = current state, parameterized bank angle of NPCG, feedback gain
% output = bank angle command

function sig_cmd = cal_QEGC(V, h, gamma, signom, KFB, auxdata)

% auxdata
S  = auxdata.S;
m  = auxdata.m;
CL = auxdata.CL;
CD = auxdata.CD;
DU = auxdata.DU;
VU = auxdata.VU;
TU = auxdata.TU;
pn = auxdata.pn;

% constants
LD = CL/CD;
beta_r = -893;

% lift
Vreal = V*VU; % m/s
hreal = h*DU; % m
rho   = cal_airdens(hreal, pn); % kg/m^3
q     = 0.5*rho*Vreal^2;
L     = q*S*CL/m; % m/s^2
L     = L/DU*TU^2; % unitless

% flight path angle of QEGC
gamma_QEGC = 1/(LD*0.5*V^2*beta_r);
% gamma_QEGC = (V^2 + L*cosd(0))/(LD*0.5*V^2*beta_r*cosd(0));

% altitude rate of QEGC
hdot_QEGC = sin(gamma_QEGC);

% current altitude rate
hdot = V*sin(gamma);

% bank anlge command compenstating QEGC
Lver = L*cos(signom) - KFB*(hdot - hdot_QEGC);
temp = Lver/L;
if temp >= 1
    temp = 1;
end
sig_cmd = acos(temp);