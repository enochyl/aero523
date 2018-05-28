function [m_slat, m_flap] = flowrate(mesh,psi)
% PURPOSE: calculates the mass flow rate between slat and main, and flap
% and main
%
% INPUTS:
%   mesh : input mesh structure, consisting of:
%     mesh.V  =  x,y coordinates of all nodes
%     mesh.E  =  triangular elements, as sequences of 3 nodes
%     mesh.IE =  interior faces and connectivities
%     mesh.BE =  boundary faces and connectivities
%   psi : stream functions
%
% OUTPUTS:
%   m_slat: mass flow rate between the slat and main
%   m_flap: mass flow rate between the flap and main
%

BE  = mesh.BE;
nBE = length(BE);

for i = 1:nBE
    if BE(i,4) == 3 % Slat
        m_slat = psi(BE(i,1));  % note: streamlines on the surface of slat or flap carry the same values.
    elseif BE(i,4) == 4 % Flap
        m_flap = -psi(BE(i,1)); % accounting for the direction with the negative sign
    end
end

