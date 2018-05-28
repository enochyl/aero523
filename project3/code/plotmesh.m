function plotmesh(mesh)
% PURPOSE: plots the provided mesh
%
% INPUTS:
%    mesh    : mesh to plot
%

figure()
t = mesh.E';
t = [t; ones(1,size(t,2))];
pdeplot(mesh.V',[],t,'mesh','on');
axis equal;
axis([-0.3 0.9 -0.35 0.25])

