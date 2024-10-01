
% ----------------------------------------------------------------------- %
% calculate derivative using the parameterized bank angle
%
%    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
%    Prof. Dae Young Lee** and Prof. Bong Wie***
%    Iowa State University, Ames, IA 50011
%    *leeyr111@gmail.com
%    **daylee@iastate.edu
%    ***bongwie@iastate.edu
% ----------------------------------------------------------------------- %
% input  = current state, guidance
% output = Jacobian

function J = cal_J(X0, sigs, e0, auxdata)

% auxdata
dsig = auxdata.dsig;

% small deviation in control variable
u1 = sigs(1) + dsig;
G1 = cal_z(X0, u1, e0, auxdata);

u2 = sigs(1) - dsig;
G2 = cal_z(X0, u2, e0, auxdata);

J = (G1(1) - G2(1))/(2*dsig);

end