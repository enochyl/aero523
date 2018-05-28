function inputcheck(mesh_size,order)
% PURPOSE: check for invalid inputs
%
% INPUTS:
%   mesh_size : size of the mesh
%   order     : order of method
%

switch mesh_size
    case 0
    case 1
    case 2
    case 3
    otherwise
        error('unsupported');
end
switch order
    case 1
    case 2
    otherwise
        error('unsupported');
end