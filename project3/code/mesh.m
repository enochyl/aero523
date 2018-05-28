function mesh = mesh(mesh_size)
% PURPOSE: This function loads the mesh file
%
% OUTPUTS:
%  mesh_size: User-defined mesh size
% OUTPUTS:
%  mesh     : the mesh
%

switch mesh_size
    case 0
        mesh.E  = load('E0.txt');
        mesh.IE = load('IE0.txt');
        mesh.V  = load('V0.txt');
        mesh.BE = load('BE0.txt');
    case 1
        mesh.E  = load('E1.txt');
        mesh.IE = load('IE1.txt');
        mesh.V  = load('V1.txt');
        mesh.BE = load('BE1.txt');
    case 2
        mesh.E  = load('E2.txt');
        mesh.IE = load('IE2.txt');
        mesh.V  = load('V2.txt');
        mesh.BE = load('BE2.txt');
    case 3
        mesh.E  = load('E3.txt');
        mesh.IE = load('IE3.txt');
        mesh.V  = load('V3.txt');
        mesh.BE = load('BE3.txt');
end


