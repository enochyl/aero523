function [uinf,aoa,gamma,c] = init
% PURPOSE: This function calculates freestream conservative state vector
%
% OUTPUTS:
%  uinf   : freestream conservative state vector
%  gamma  : specific heat ratio
%  aoa    : angle of attack
%  c      : chord length
%

Minf  = 0.25;            % Freestream Mach number
gamma = 1.4;             % Specific heat ratio
aoa   = 5*pi/180;        % Angle of attack
c = 0.5588;              % Chord length

uinf(1,1) = 1;
uinf(1,2) = Minf*cos(aoa);
uinf(1,3) = Minf*sin(aoa);
uinf(1,4) = 1/(gamma-1)/gamma + Minf^2/2;

