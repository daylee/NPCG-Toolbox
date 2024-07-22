
function J = cal_J_G(X0, convar, e0, auxdata)

% auxdata
BAS  = auxdata.BAS;
dsig = auxdata.dsig;

if BAS <= 4

    % small deviation in control variable
    u1 = convar(1) + dsig;
    G1 = cal_G(X0, u1, e0, auxdata);

    u2 = convar(1) - dsig;
    G2 = cal_G(X0, u2, e0, auxdata);

    J = (G1(1) - G2(1))/(2*dsig);

else

    % small deviation in control variable
    u1 = [convar(1) + dsig, convar(2)];
    G1 = cal_G(X0, u1, e0, auxdata);

    u2 = [convar(1) - dsig, convar(2)];
    G2 = cal_G(X0, u2, e0, auxdata);

    u3 = [convar(1), convar(2) + dsig];
    G3 = cal_G(X0, u3, e0, auxdata);

    u4 = [convar(1), convar(2) - dsig];
    G4 = cal_G(X0, u4, e0, auxdata);

    % Jacobian components
    j11 = (G1(1) - G2(1))/(2*dsig);
    j12 = (G3(1) - G4(1))/(2*dsig);

    % Jacobian, 1 x 2
    J = [j11 j12];

end