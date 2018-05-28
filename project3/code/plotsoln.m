function plotsoln(mesh, un, sstring, order)
% PURPOSE: Plots cell data U, post-processed to compute scalar sname,
% over a given mesh
%
% INPUTS:
%    mesh    : mesh on which the solution was computed
%    un      : solution
%    sstring : state to be plotted
%    order   : order of the method used
%

s = zeros(size(un,1),1);
gamma = 1.4;
rho  = un(:,1);
rhou = un(:,2);
rhov = un(:,3);
rhoE = un(:,4);

switch lower(sstring)
    case 'mach'
        v = sqrt(rhou.^2 + rhov.^2)./rho;
        p = (gamma - 1)*(rhoE - 0.5*rho.*v.^2);
        c = sqrt(gamma*p./rho);
        s = v./c;
    case 'pressure'
        v = sqrt(rhou.^2 + rhov.^2)./rho;
        s = (gamma - 1)*(rhoE - 0.5*rho.*v.^2);
    otherwise
        error 'unsupported';
end

figure()
colormap(jet(2^5));
t = mesh.E';
t = [t; ones(1,size(t,2))];
pdeplot(mesh.V', [],t,'xydata',s,'xystyle','interp', 'colormap', 'jet'); 
hold on;
% pdeplot(mesh.V',[],t,'mesh','on');
hold off
axis equal; axis([-0.3 0.9 -0.35 0.25]); axis off;
switch lower(sstring)
    case 'mach'
        switch order
            case 1
                caxis([0 0.6]);
                shading flat
            case 2
                caxis([0 0.6]);
        end
    case 'pressure'
        switch order
            case 1
                caxis([0.62 0.76]);
                shading flat
            case 2
                caxis([0.62 0.76]);
        end
end
switch order
    case 1
        shading flat
    case 2
end


