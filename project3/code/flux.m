function [F, smag] = flux(UL, UR, n)
% PURPOSE: This function calculates the flux for the Euler equations
% using the Roe flux function
%
% INPUTS:
%    UL: conservative state vector in left cell
%    UR: conservative state vector in right cell
%     n: normal pointing from the left cell to the right cell
%
% OUTPUTS:
%  F   : the flux out of the left cell (into the right cell)
%  smag: the maximum propagation speed of disturbance
%

gamma = 1.4;
gmi = gamma-1.0;

% process left state
rL = UL(1);
uL = UL(2)/rL;
vL = UL(3)/rL;
unL = uL*n(1) + vL*n(2);
qL = sqrt(UL(2)^2 + UL(3)^2)/rL;
pL = (gamma-1)*(UL(4) - 0.5*rL*qL^2);
if ((pL<=0) || (rL<=0)), error 'Non-physical state!', end;
rHL = UL(4) + pL;
HL = rHL/rL;
cL = sqrt(gamma*pL/rL);

% left flux
FL = UL; % for allocation
FL(1) = rL*unL;
FL(2) = UL(2)*unL + pL*n(1);
FL(3) = UL(3)*unL + pL*n(2);
FL(4) = rHL*unL;

% process right state
rR = UR(1);
uR = UR(2)/rR;
vR = UR(3)/rR;
unR = uR*n(1) + vR*n(2);
qR = sqrt(UR(2)^2 + UR(3)^2)/rR;
pR = (gamma-1)*(UR(4) - 0.5*rR*qR^2);
if ((pR<=0) || (rR<=0)), error 'Non-physical state!', end;
rHR = UR(4) + pR;
HR = rHR/rR;
cR = sqrt(gamma*pR/rR);

% right flux
FR = UR; % for allocation
FR(1) = rR*unR;
FR(2) = UR(2)*unR + pR*n(1);
FR(3) = UR(3)*unR + pR*n(2);
FR(4) = rHR*unR;

% difference in states
du = UR - UL;

% Roe average
di     = sqrt(rR/rL);
d1     = 1.0/(1.0+di);

ui     = (di*uR + uL)*d1;
vi     = (di*vR + vL)*d1;
Hi     = (di*HR + HL)*d1;

af     = 0.5*(ui*ui+vi*vi );
ucp    = ui*n(1) + vi*n(2);
c2     = gmi*(Hi - af);
if (c2<=0), error 'Non-physical state!', end;
ci     = sqrt(c2);
ci1    = 1.0/ci;

% eigenvalues
z = zeros(3,1);
l(1) = ucp+ci;
l(2) = ucp-ci;
l(3) = ucp;

% entropy fix
epsilon = ci*.1;
for i=1:3,
  if ((l(i)<epsilon) && (l(i)>-epsilon)),
    l(i) = 0.5*(epsilon + l(i)*l(i)/epsilon);
  end
end

l = abs(l); l3 = l(3);

% average and half-difference of 1st and 2nd eigs
s1    = 0.5*(l(1) + l(2));
s2    = 0.5*(l(1) - l(2));

% left eigenvector product generators (see Theory guide)
G1    = gmi*(af*du(1) - ui*du(2) - vi*du(3) + du(4));
G2    = -ucp*du(1)+du(2)*n(1)+du(3)*n(2);

% required functions of G1 and G2 (again, see Theory guide)
C1    = G1*(s1-l3)*ci1*ci1 + G2*s2*ci1;
C2    = G1*s2*ci1          + G2*(s1-l3);
 
% flux assembly
F = FL; % for allocation
F(1)    = 0.5*(FL(1)+FR(1))-0.5*(l3*du(1) + C1   );
F(2)    = 0.5*(FL(2)+FR(2))-0.5*(l3*du(2) + C1*ui + C2*n(1));
F(3)    = 0.5*(FL(3)+FR(3))-0.5*(l3*du(3) + C1*vi + C2*n(2));
F(4)    = 0.5*(FL(4)+FR(4))-0.5*(l3*du(4) + C1*Hi + C2*ucp  );

% max wave speed
smag = max(l);
