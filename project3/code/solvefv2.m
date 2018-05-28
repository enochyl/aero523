function [un,iter,Rmax,cl,cd,cl_conv,cd_conv,xpos,cpd] = solvefv2(mesh)
% PURPOSE: solves the compressible, subsonic flow over a high-lift airfoil
% using second-order finite volume method
%
% INPUTS:
%   mesh : input mesh structure
%     mesh.V  =  x,y coordinates of all nodes
%     mesh.E  =  triangular elements, as sequences of 3 nodes
%     mesh.IE =  interior faces and connectivities
%     mesh.BE =  boundary faces and connectivities
%
% OUTPUTS:
%   un          : solution
%   iter        : iterations
%   Rmax        : residual convergence
%   cl          : component lift coefficients
%     cl.main = main lift coefficient
%     cl.slat = slat lift coefficient
%     cl.flap = flap lift coefficient
%   cd          : component drag coefficients
%     cd.main = main drag coefficient
%     cd.slat = slat drag coefficient
%     cd.flap = flap drag coefficient
%   cl_conv     : lift coefficient convergence vector
%   cd_conv     : drag coefficient convergence vector
%   xpos        : x-position of boundary edge midpoints
%   cpd         : surface pressure coefficient distribution
%

%% Constants definitions
Rtol = -8; % Residual tolerance
Rmax = 0;
CFL = 0.9;
iter = 0;

%% Extracting quantities from mesh
V = mesh.V;  
E = mesh.E;
IE = mesh.IE;
BE = mesh.BE;

%% Number of elements, interior edges, and boundary edges

nIE   = length(IE);
nBE   = length(BE);

%% Interior edge normal and length
% Coordinates of interior nodes
ix1 = V(IE(1:nIE,1),1); ix2 = V(IE(1:nIE,2),1); 
iy1 = V(IE(1:nIE,1),2); iy2 = V(IE(1:nIE,2),2);

% Edge length
s.IElength = sqrt((ix1-ix2).^2 + (iy1-iy2).^2);

% Normal components
s.IEnormal(:,1) =  (iy2-iy1)./s.IElength;
s.IEnormal(:,2) = -(ix2-ix1)./s.IElength;

%% Boundary edge normal and length
% Coordinates of boundary nodes
bx1 = V(BE(1:nBE,1),1); bx2 = V(BE(1:nBE,2),1); 
by1 = V(BE(1:nBE,1),2); by2 = V(BE(1:nBE,2),2);

% Edge length
s.BElength = sqrt((bx1-bx2).^2 + (by1-by2).^2);

% Normal components
s.BEnormal(:,1) =  (by2-by1)./s.BElength;
s.BEnormal(:,2) = -(bx2-bx1)./s.BElength;

%% Cell area %%
% A = 1/2(x1y2 + x2y3 + x3y1 - x1y3 - x2y1 - x3y2)
area  = 1/2*abs((V(E(:,1),1).*V(E(:,2),2) + V(E(:,2),1).*V(E(:,3),2)...
        + V(E(:,3),1).*V(E(:,1),2) - V(E(:,1),1).*V(E(:,3),2)...
        - V(E(:,2),1).*V(E(:,1),2) - V(E(:,3),1).*V(E(:,2),2)));
geo.area = repmat(area,1,4);

%% Centroids %%
cx = 1/3*(V(E(:,1),1)+V(E(:,2),1)+V(E(:,3),1)); % Centroid x-coordinates
cy = 1/3*(V(E(:,1),2)+V(E(:,2),2)+V(E(:,3),2)); % Centroid y-coordinates
geo.cent = [cx,cy];

%% Initial guess and/or Freestream conditions
[uinf,aoa,gamma,c] = init;  % Initialize u0
rhoinf = uinf(1);   % Freestream density
uxinf = uinf(2)/rhoinf; % Freestream x-velocity
uyinf = uinf(3)/rhoinf; % Freestream y-velocity
rhoEinf = uinf(4);  % Energy
Uinf = sqrt(uxinf^2+uyinf^2);   % Freestream velocity
Pinf = (gamma-1)*(rhoEinf - 1/2*rhoinf*Uinf^2); % Freestream pressure
q = 1/2*rhoinf*Uinf^2;  % Free stream dynamic pressure

[un] = solvefv1(mesh());

%% Second-order finite-volume
while Rmax > Rtol
    
    [R,smag] = residualfv2(mesh,un,s,uinf,geo);

    % Time step
    ldt = 2.*CFL./smag;
    v_ldt = repmat(ldt,1,4);
    
    %% RK2 Scheme
    % Intermediate state
    unFE = un - R.*v_ldt;
    [RFE,~] = residualfv2(mesh,unFE,s,uinf,geo);
    % State update
    unp1 = 1/2*(un + unFE - RFE.*v_ldt);
    un = unp1;
    
    iter = iter + 1;
    Rmax(iter) = log10((max(max(abs(RFE)))));
    
    if mod(iter,10)==0
        fprintf(1, 'n = %d, log(Rmax) = %.5f,\n', iter, Rmax(iter));
    end
    
    %% Cl and Cd calculations at each iteration
    cl.main = 0; cl.slat = 0; cl.flap = 0;
    cd.main = 0; cd.slat = 0; cd.flap = 0;

    for i = 1:nBE
        if BE(i,4) ~= 1 % Ignoring freestream boundaries
            eL = BE(i,3); UL = un(eL,:); 
            % Boundary Pressure
            [~, ~, P] = wallflux(UL,s.BEnormal(i,:));
            % Forces
            Fx = P*s.BElength(i)*s.BEnormal(i,1); Fy = P*s.BElength(i)*s.BEnormal(i,2); 
            % Lift and Drag
            l = Fy*cos(aoa) - Fx*sin(aoa); d = Fy*sin(aoa) + Fx*cos(aoa);
            switch BE(i,4)
                case 2
                    cl.main = cl.main + l/q/c;
                    cd.main = cd.main + d/q/c;
                case 3
                    cl.slat = cl.slat + l/q/c;
                    cd.slat = cd.slat + d/q/c;
                case 4
                    cl.flap = cl.flap + l/q/c;
                    cd.flap = cd.flap + d/q/c;
            end
        end
    end
    cl_conv(iter) = cl.main + cl.slat + cl.flap;
    cd_conv(iter) = cd.main + cd.slat + cd.flap;
end
%% Cp distribution
j = 1;
for i = 1:nBE
    if BE(i,4) ~= 1 % Ignore freestream boundaries
        eL = BE(i,3); UL = un(eL,:); 
        [~, ~, pL] = wallflux(UL,s.BEnormal(i,:));
        cpd(j) = (pL-Pinf)/q;
        xpos(j) = (V(BE(i,1),1)+V(BE(i,2),1))/2;
        j = j + 1;
    end
end

