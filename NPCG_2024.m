% NPCG_2024.m
% To be released  in September 2024 for AIAA Online Course
% Numerical Predictor Corrector Guidance Toolbox (Version 2024)... 
% with various options of a standard linear bank angle function...
% and other variants of bank-angle shaping functions, as applied to...
% a variety of atmospheric entry guidance problems

% Coded by Dr. Youngro Lee for his PhD work under the supervision of...  
% Prof. Dae Young Lee (daylee@iastate.edu) and Prof. Bong Wie (bongwie@iastate.edu) 
% Iowa State University, Ames, IA 50011


clear;clc;close all

% add a path containing commonly used functions
addpath('comfncs')

% ----------------------------------------------------------------------- %
% select the entry example case / define design parameters
% ----------------------------------------------------------------------- %
%  Case_Number = atmospheric entry examples
%     1: Mars robotic mission (Mars Science Laboratory model)
%     2: Mars human mission model of NASA (2018) 
%     3: Mars human mission model updated in  Signialo et al (2024)
%     4: Apollo 10
%     5: CAV-H
%
% GAT = guidance activation time in seconds since entry interface (EI)
%
% BAL = bank angle limit
%     1: simple magnitude contraint
%     2: rate and acceleration constraints
%     siglmt   = angle limit
%     sigdlmt  = rate limit
%     sigddlmt = acceleration limit
%
% BAS = bank angle shaping options
%     1: linear function (Lu 2014)
%     2: exponential function (Liang and Zhu 2021)
%     3: exponential function (Youngro Lee)
%     4: logistic function (Youngro Lee)
%     5: two-segment linear function (Youngro Lee)  
%     6: quadratic function(Li et al., 2019)
%     7: automated linear function (Youngro Lee)
%  Check lines 195-227 for additional parameters to be adjusted for each BAS option  

% BRL = bank reversal logic
%     1: conventional logic by Tu, Munir, Mease, and Bayard (JGCD 2000) 
%     2: predictive logic by K.M. Smith (2016)
%     KBR = damping ratio of the predictive lateral guidance
% ----------------------------------------------------------------------- %

Case_Number = 3;


 
switch   Case_Number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1

        Mars_robotic

       GAT= 60;

        BAL = 1;
        siglmt = 180*d2r;   % bank angle limit
        sigdlmt =  15*d2r;  % bank rate limit
        sigddlmt = 15*d2r;  % bank acc limit

        sig0=70*d2r; sigf =20*d2r;

        BAS = 1;   % 7.98 km
        BA= 2; % 7.98 km
        BAS= 3; % 7.78 km
       % BAS = 4;   % 8.05 km
       % BAS =5; % 8.18 km
        BAS = 6; % 8.07 km
        %BAS = 7;   % 8.15 km
       
        BRL = 2; KBR = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 

         Mars_human
         GAT = 160 ;
        BAL = 1;
        siglmt = 180*d2r;
        sigdlmt =  15*d2r;  % bank rate limit
        sigddlmt = 15*d2r;  % bank acc limit
   
        sig0=70*d2r; sigf =20*d2r; 

  BAS = 1;    % crashed
  % BAS= 2; % 0.77 km
 % BAS= 3; %1.5 km
 % BAS = 4;     %1.45 km
 % BAS = 5;  % 2.72 km
 % BAS= 6; % 3.45 km
 % BAS = 7;    %1.9 km

      BRL = 2; KBR = 7;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 3

         Mars_human_new
    
        BAL = 1;
        siglmt = 180*d2r;    % bank angle limit
        sigdlmt =   20*d2r;  % bank rate limit
        sigddlmt = 20*d2r;  % bank acc limit
  
        sig0=70*d2r; sigf =40*d2r;

  
  BAS = 1;  GAT=175;    % 1.11 km  
  
  BRL = 2; KBR = 7;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 4

         Earth_Apollo10

        GAT = 90;

         BAL = 1;
        siglmt = 180*d2r;
        sigdlmt =  15*d2r;  % bank rate limit
        sigddlmt = 15*d2r;  % bank acc limit

          sig0=70*d2r; sigf =20*d2r;

         BAS= 1;   % 5.99 km  
        % BAS = 2; % failed
        % BAS = 3; % failed
        % BAS = 4;   % failed
       % BAS = 5;  % 5.65 km
       % BAS = 6; % failed
       % BAS = 7; % 5.68 km

           BRL = 2; KBR = 10;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5

         Earth_CAVH

        GAT = 300;

        BAL = 1;
        siglmt = 180*d2r;
        sigdlmt =  15*d2r;  % bank rate limit
        sigddlmt = 15*d2r;  % bank acc limit
        
          sig0=70*d2r; sigf =20*d2r;

        BAS = 1; % failed
        BAS = 2; % 27.58 km
        BAS = 3; % failed
        BAS = 4; % 28.82 km
        % BAS = 5; % not applicable
        % BAS = 6; % not applicable
        % BAS =7; % not applicable

           BRL = 2;  KBR = 5;

        POD = 1; % 0 = no damping; 1 = active damping of phugoid oscillations 
        VFB = 4500; % velocity threshold
        KFB = 15; % feedback gain of active damping of phugoid oscillations
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bank angle shaping (BAS) options

if BAS == 1 % linear function (Lu 2014)    
    sigs = sig0;    
    auxdata.sigf = sigf;

elseif BAS == 2 % exponential function (Liang and Zhu 2021)
    KEF  = 1;      % design parameter
    sigs = sig0;       
    auxdata.sigf = sigf;    
    auxdata.KEF = KEF;

elseif BAS == 3 % exponential function (Youngro Lee)
    KEF = 0.6;     % design parameter 
    sigs = sig0;
    auxdata.KEF = KEF;

elseif BAS == 4 % logistic function (Youngro Lee)
    KLF = 1.3;  % design parameter     
    sigs = sig0;
    auxdata.KLF = KLF;

elseif BAS == 5 % two-segment linear function  
    sigm = 60*d2r; % initial guess of parameter to be found
    sigs = [sig0; sigm];
    auxdata.sigf = sigf;

elseif BAS == 6 % quadratic function(Li et al., 2019)
    sigm = 100*d2r; % initial guess of parameter to be found
    sigs = [sig0; sigm];
    auxdata.sigf = sigf;

elseif BAS == 7 % automated linear function (Youngro Lee)
    sigs = [sig0; sigf];
end

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
    BRprev =   - 1;
else
    BRprev =    1;
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
auxdata.BAS      = BAS;
auxdata.KBR      = KBR;
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
        bankcmd =0;

    else % guidance begin

        inum = 0;
        % numerical predictor-corrector algorithm execution
        while 1

            % range calculation
            G = cal_G(X0, sigs, ec, auxdata);

            % range condition
            if abs(G) <= 1e2/DU
                break
            end

            % Jacobian
            J = cal_J_G(X0, sigs, ec, auxdata);

            % psuedo inverse
            invJ = pinv(J);

            % Newton–Raphson method
            sigs = sigs - G*invJ; % solution update

            inum = inum + 1;
            fprintf('Correction at t = %3.0f s, h = %3.0f km , V = %3.0f m/s \n',...
                tc*TU, hc*DU*1e-3, Vc*VU )

            % numerical condition
            if norm(J) < 1e-6
                break
            end
        end
        testdata(end+1,:) = [tc ec sigs' inum];

        if Case_Number == 5 && POD == 1
            if Vc >= VFB/VU
                sigs = cal_QEGC(Vc, hc, gamc, sigs, KFB, auxdata);
            end
        end

        if BAL == 1
            % bank angle  limit
            if abs(sigs(1)) >= siglmt
                sigs(1) = siglmt;
            end
        end

        % bank reversal logic
        if BRL == 1
            BR = Aux_BR(X0, BRprev, ex_case, auxdata);
        elseif BRL == 2
            BR = Aux_BR_prdt(X0, ec, sigs, BRprev, auxdata);
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
    [tt,XX] = ode45(@(t,X) Aux_dyn(t, X, bankcmd, auxdata), tspan, X0);

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

% energy
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
disp(['Targeting Error: ', num2str(Rgof),' km'])
disp(['Simulation Time: ', num2str(elpdtime), ' s'])
disp(' ')

% plotting simulation variables
%plotting
myplot
