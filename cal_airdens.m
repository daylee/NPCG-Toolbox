
% ----------------------------------------------------------------------- %
% calculate air denstiy
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %
% input  = altitude in m
% output = air density in kg/m^3

function rho = cal_airdens(h, pn)

switch pn % planet, either Earth or Mars
    case 1
        % Earth atmosphere (hicks 2009)
        beta = 0.14; % 1/km
        rhos = 1.225;
        rho = rhos*exp(-beta*h*1e-3);

    case 2
        % Mars atmosphere (AAS 18-485)
        beta = -0.000105 * h;
        beta1 = 559.351005946503;
        beta2 = 188.95110711075;

        T = (1.4e-13)*(h.^3) - (8.85e-9)*(h.^2) - (1.245e-3)*h + 205.3645;

        rho0 = beta1./(beta2.*T);
        rho  = rho0.* exp(beta);
end