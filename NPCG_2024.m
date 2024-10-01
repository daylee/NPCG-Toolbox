
% ----------------------------------------------------------------------- %
%      Numerical Predictor-Corrector Guidance Toolbox (Version 2024)
% ----------------------------------------------------------------------- %
%    To be released in September 2024 for AIAA Online Course
%    Numerical Predictor-Corrector Guidance Toolbox...
%    with a standard linear bank angle parameterization...
%    and other variants of bank angle parameterizations, as applied to...
%    a variety of atmospheric entry guidance problems
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %
%
%
% ----------------------------------------------------------------------- %
%  select the entry example (Case_Number) out of five cases
%  select the bank angle parameterization (BAP) out of seven options
% ----------------------------------------------------------------------- %
% Case_Number = atmospheric entry examples
%     1: Mars robotic mission (Mars Science Laboratory model)
%     2: Mars manned mission model in (Jiang et al., 2019)
%     3: Mars manned mission model updated in (Signialo et al, 2024)
%     4: Apollo 10 (Szelc 1969)
%     5: CAV-H (Lu 2013)
%
% BAP = bank angle parameterization
%     1: linear function            (Lu 2014)
%     2: exponential function       (Liang and Zhu 2021)
%     3: exponential function       (Youngro Lee)
%     4: logistic function          (Youngro Lee)
%
% Note!!! each BAP option has guidance parameters to be adjusted, which
%         can be found in "ex_" scripts.
%
% ----------------------------------------------------------------------- %
%  GAT, BRL, and BAL should be wisely selected for a good guidance
%  performance. A good parameter set for each combination are already given
%  in "ex_***", and they can be adjusted according to the user's desire.
% ----------------------------------------------------------------------- %
%
% GAT = guidance activation time in seconds since entry interface (EI)
%
% BRL = bank reversal logic
%     1: conventional logic by Tu, Munir, Mease, and Bayard (JGCD 2000)
%           dlpsT = heading angle error thresholds
%     2: predictive logic by K.M. Smith (AAS 2016)
%           KBR = damping ratio of the predictive lateral guidance
%
% BAL = bank angle constraints
%     1: simple magnitude contraint
%           siglmt   = magnitude limit
%     2: rate and acceleration constraints
%           sigdlmt  = rate limit
%           sigddlmt = acceleration limit
%
% ----------------------------------------------------------------------- %

clear;clc;close all

% entry example cases
Case_Number = 1;

% bank angle parameterization options
BAP = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch Case_Number

    case 1

        ex_Mars_robotic

    case 2

        ex_Mars_manned

    case 3

        ex_Mars_manned_new

    case 4

        ex_Earth_Apollo10

    case 5

        ex_Earth_CAVH

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalization
DU = rp;
VU = sqrt(mu/DU);
TU = DU/VU;
rp = rp/DU;
mu = mu*DU^-3*TU^2;
w  = w*TU;
Vf = Vf/VU;
rf = rf/DU;
V0 = V0/VU;
r0 = r0/DU;

% initial state vector
X0 = [r0, theta0, phi0, V0, gamma0, psi0];

% initial and final energy conditions
e0 = mu/r0 - (V0^2)/2;
ef =  mu/rf - (Vf^2)/2;

% initial bank reversal value
delpsi0 = cal_delpsi(X0, thetatgt, phitgt);
if sign(delpsi0) > 0
    BRprev = -1;
else
    BRprev = 1;
end

% increment for Newton–Raphson method
dsig = 0.01*d2r;

% propagation time step and guidance frequency
tstp = 1; % s

% ode45 setting
opt = odeset('RelTol',1e-6,'AbsTol',1e-9);

% auxdata
auxdata.pn       = pn;
auxdata.rp       = rp;
auxdata.mu       = mu;
auxdata.w        = w;
auxdata.S        = S;
auxdata.m        = m;
auxdata.CL       = CL;
auxdata.CD       = CD;
auxdata.DU       = DU;
auxdata.VU       = VU;
auxdata.TU       = TU;
auxdata.theta0   = theta0;
auxdata.phi0     = phi0;
auxdata.thetatgt = thetatgt;
auxdata.phitgt   = phitgt;
auxdata.ef       = ef;
auxdata.rf       = rf;
auxdata.BAP      = BAP;
auxdata.KBR      = KBR;
auxdata.dlpsT    = dlpsT;
auxdata.opt      = opt;
auxdata.dsig     = dsig;

% initialization
tc = 0; % current time
ec = e0; % current energy
hc = r0 - rp; % current altitude
Vc = V0; % current velocity
sigcmdprv = 0; % previous step bank angle command
sigcmdprv2 = 0; % two-step previous bank angle command

% data array
X = X0;
T = tc;
SIG = 0;
testdata = [];

% ------------------- software-in-the-loop Simulation ------------------- %
tic
while ec < ef

    if tc <= GAT/TU
        bankcmd = 0;

    else % guidance begin

        inum = 0;
        % numerical predictor-corrector algorithm execution
        while 1

            % range calculation
            G = cal_z(X0, sigs, ec, auxdata);

            % range condition
            if abs(G) <= 1e3/DU
                break
            end

            % derivative
            J = cal_J(X0, sigs, ec, auxdata);

            % Newton–Raphson method
            sigs = sigs - G/J; % solution update

            inum = inum + 1;
            fprintf('   Correction at t = %3.0f s \n', tc*TU)

            % numerical condition
            if norm(J) < 1e-6
                break
            end
        end
        testdata(end+1,:) = [tc ec sigs' inum];

        if BAL == 1
            % bank angle  limit
            if abs(sigs(1)) >= siglmt
                sigs(1) = siglmt;
            end
        end

        % bank reversal logic
        if BRL == 1
            BR = cal_BR(X0, BRprev, Case_Number, auxdata);
        elseif BRL == 2
            BR = cal_BR_prdt(X0, ec, sigs, BRprev, auxdata);
        end

        if BR ~= BRprev
            Rgo = DU*1e-3*cal_sphdis(X0(2), X0(3), thetatgt, phitgt);
            fprintf('       Bank reversal at Rgo = %3.0f km \n', Rgo )
        end
        BRprev = BR;

        % bank command
        bankcmd = sigs(1)*BR;

        if BAL == 2
            % bank angle magnitude limit
            % calculation of lower bound

            tmp1 = sigcmdprv - sigdlmt*tstp; % from rate
            tmp2 = 2*sigcmdprv - sigcmdprv2 - sigddlmt*tstp^2; % from rate
            minbnk = max([tmp1, tmp2]);

            % calculation of upper bound
            tmp3 = sigcmdprv + sigdlmt*tstp; % from rate
            tmp4 = 2*sigcmdprv - sigcmdprv2 + sigddlmt*tstp^2; % from rate
            maxbnk = min([tmp3, tmp4]);

            if bankcmd >= maxbnk
                bankcmd = maxbnk;
            elseif bankcmd <= minbnk
                bankcmd = minbnk;
            end
            % bank angle limit
            if abs(bankcmd) >= siglmt
                bankcmd = siglmt;
            end
            sigcmdprv2 = sigcmdprv;
            sigcmdprv = bankcmd;
        end

    end

    tspan = [tc, tc + tstp/TU];
    [tt,XX] = ode45(@(t,X) aux_dyn(t, X, bankcmd, auxdata), tspan, X0);

    Xn = XX(end,:);
    tc = tt(end);

    % for the next step
    hc = Xn(1) - rp;
    Vc = Xn(4);
    gamc = Xn(5);
    X0 = Xn;
    ec = mu/Xn(1) - (Xn(4)^2)/2;

    % save data
    X(end+1,:) = Xn;
    T(end+1,1) = tc;
    SIG(end+1,1) = bankcmd;

    % guidance termination
    Rgo = DU*cal_sphdis(X0(2), X0(3), thetatgt, phitgt);
    if Rgo < 100 % m
        disp('  Rgo < 100 m')
        break
    end

end
elpdtime = toc;

% --------------------------- post-processing --------------------------- %
% state variables
r = X(:,1); theta = X(:,2); phi = X(:,3);
V = X(:,4); gamma = X(:,5); psi = X(:,6);

% energy-like variable
e =  mu./r - 0.5*V.^2;

% time
time = T*TU;

% bank angle
sigcmd = SIG;
sigmadata = [time, sigcmd];

% downrange and crossrange
nn = length(time); DR = zeros(nn,1); CR = DR; delpsi = DR;
for ii = 1:nn
    [dr, cr] = cal_drdc(theta0, phi0, thetatgt, phitgt, theta(ii), phi(ii));
    DR(ii) = dr;
    CR(ii) = cr;
    delpsi(ii) = cal_delpsi(X(ii,:), thetatgt, phitgt);
end
DR0 = cal_sphdis(theta0, phi0, thetatgt, phitgt); % mission range-to-go
DRtogo = (DR0 - DR)*DU; % downrange-to-go
Rgof = cal_sphdis(theta(end), phi(end), thetatgt, phitgt)*DU*1e-3;

% path constraints
Vreal = V*VU; % m/s
hreal = (r - rp)*DU; % m
rho = cal_airdens(hreal, pn); % kg/m^3
q = 0.5*rho.*Vreal.^2; % Pa
D = q.*CD.*S/m; % m/s^2
L = q.*CL.*S/m; % m/s^2
A = sqrt(L.^2 + D.^2)/9.81; % g
qdot = kq*(rho.^kqN).*(Vreal.^kqM); % kW/m^2

% unit conversion
h      = (r-rp)*DU*1e-3;
V      = V*VU*1e-3;
theta  = theta*r2d;
phi    = phi*r2d;
gamma  = gamma*r2d;
psi    = psi*r2d;
DR     = DR*DU*1e-3;
CR     = CR*DU*1e-3;
DRtogo = DRtogo*1e-3;

% fianl targeting error and simulation time
disp(' ')
disp([' Simulation Time: ', num2str(elpdtime, 3), ' sec'])
disp(' ')

% plotting simulation variables
myplot