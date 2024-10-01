
% ----------------------------------------------------------------------- %
% calculate downrange and crossrange
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %
% input  = longitude and latitude of initial, target, and current
%          locations in rad
% output = downrange and crossrange in rad

function [DR,CR] = cal_drdc(theta0, phi0, thetaf, phif, theta, phi)

% sin/cos of angles
cth0 = cos(theta0); sth0 = sin(theta0); cph0 = cos(phi0); sph0 = sin(phi0);
cthf = cos(thetaf); sthf = sin(thetaf); cphf = cos(phif); sphf = sin(phif);
cth  = cos(theta);  sth  = sin(theta);  cph  = cos(phi);  sph  = sin(phi);

% unit distance from X0 to X on unit sphere
c = cal_sphdis(theta0, phi0, theta, phi);

% X0 = initial position vector
% Xt = target position vector
% X  = current position vector on the unit sphere
X0 = [ cth0*cph0;  sth0*cph0; sph0];
Xt = [ cthf*cphf;  sthf*cphf; sphf];
X  = [ cth*cph;  sth*cph; sph];

% normal vector of X0 and Xt
nX0Xt = cross(X0, Xt);
nX0Xt = nX0Xt/norm(nX0Xt);

% normal vector of X0 and X
if theta0 == theta && phi0 == phi % initial point
    nX0X = nX0Xt;
else
    nX0X = cross(X0, X);
    nX0X = nX0X/norm(nX0X);
end

% included angle between X0toX and X0toXt
A = acos(dot(nX0X,nX0Xt));
A = real(A);

% crossrange; sin(a) = sin(A)sin(c)
a = asin(sin(A)*sin(c));
CR = a;

% downrange; cos(c) = cos(a)cos(b)
b = acos(cos(c)/cos(a));
DR = b;

temp0 = cross(nX0X, nX0Xt);
temp1 = temp0/norm(temp0);

% define crossrange sign
temp2 = dot(temp1,X0);

if sign(temp2) < 0
    CR = -CR;
end