
% ----------------------------------------------------------------------- %
% calculate a great-circle distance on the unit spherical surface
% ----------------------------------------------------------------------- %
% input  = longitude and latitude of two locations in rad
% output = angle distance in rad

function ang = cal_sphdis(theta1, phi1, theta2, phi2)

ang = acos(sin(phi2)*sin(phi1) + cos(phi2)*cos(phi1)*cos(theta2-theta1));