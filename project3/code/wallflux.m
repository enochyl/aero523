function [FL, smag, pL] = wallflux(UL, n)
% PURPOSE: This function calculates the flux for the Euler equations
% at the wall as well as wall pressure.
%
% INPUTS:
%    UL: conservative state vector in left cell
%     n: normal pointing from the left cell to the right cell
%
% OUTPUTS:
%  F   : the flux out of the left cell (into the right cell)
%  smag: the maximum propagation speed of disturbance
%  pL  : pressure at the wall
%

gamma = 1.4;

rL  = UL(1);
uL  = UL(2)/rL;
vL  = UL(3)/rL;
unL = uL*n(1) + vL*n(2);
qL  = sqrt(UL(2)^2 + UL(3)^2)/rL;
utL = sqrt(qL^2 - unL^2);
pL = (gamma-1)*(UL(4) - 0.5*rL*utL^2);
if ((pL<=0) || (rL<=0)), error 'Non-physical state!', end;
cL = sqrt(gamma*pL/rL);

smag = abs(unL) + cL;

FL = zeros(size(UL));
FL(2) = pL*n(1);
FL(3) = pL*n(2);