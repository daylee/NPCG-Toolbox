
% ----------------------------------------------------------------------- %
% calculate a great-circle distance on an unit spherical surface
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %
% input  = longitude and latitude of two locations in rad
% output = angle distance in rad

function ang = cal_sphdis(theta1, phi1, theta2, phi2)

ang = acos(sin(phi2)*sin(phi1) + cos(phi2)*cos(phi1)*cos(theta2-theta1));