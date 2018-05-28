function [R,smag] = residualfv2(mesh,un,s,uinf,geo)
% PURPOSE: computes the residual using second-order fluxes
%
% INPUTS:
%   mesh : input mesh structure, consisting of:
%     mesh.V  =  x,y coordinates of all nodes
%     mesh.E  =  triangular elements, as sequences of 3 nodes
%     mesh.IE =  interior faces and connectivities
%     mesh.BE =  boundary faces and connectivities
%   un   : input state
%   s    : edge normals and lengths
%     s.IEnormal = interior edge normal
%     s.IElength = interior edge length
%     s.BEnormal = boundary edge normal
%     s.BElength = boundary edge length
%   unif : freestream states
%   geo  : element gemoetries
%     geo.area = cell area
%     geo.cent = cell centroid
%
% OUTPUTS:
%   R          : residual
%   smag       : the maximum propagation speed of disturbance
%

%% Extracting quantities from inputs
V = mesh.V;
E = mesh.E;
IE = mesh.IE;
BE = mesh.BE;

nelem = length( E);
nIE   = length(IE);
nBE   = length(BE);

IEnormal = s.IEnormal;    % Interior edge normal
IElength = s.IElength;    % Interior edge length
BEnormal = s.BEnormal;    % Boundary edge normal
BElength = s.BElength;    % Boundary edge length

area = geo.area;    % Cell area
cent = geo.cent;    % Centroid

%% Residual calculation (2nd order)
gradx = zeros(nelem,4); % Initialization of the gradient, x-comp
grady = zeros(nelem,4); % Initialization of the gradient, y-comp

for i = 1:nIE
    eL = IE(i,3); eR = IE(i,4);

    % Gradient Calculations (interior)
    uhat = 1/2.*(un(eL,:) + un(eR,:));    
    grad_comp = uhat'*IEnormal(i,:).*IElength(i);

    gradx(eL,:) = gradx(eL,:) + grad_comp(:,1)'; grady(eL,:) = grady(eL,:) + grad_comp(:,2)';   % Left cell
    gradx(eR,:) = gradx(eR,:) - grad_comp(:,1)'; grady(eR,:) = grady(eR,:) - grad_comp(:,2)';   % Right cell
end
for i = 1:nBE
    eL = BE(i,3);

    % Gradient Calculations (boundary)
    if BE(i,4) == 1;    % Freestream edges
        uhat = 1/2*(un(eL,:) + uinf);
        grad_comp = uhat'*BEnormal(i,:).*BElength(i);
        gradx(eL,:) = gradx(eL,:) + grad_comp(:,1)'; 
        grady(eL,:) = grady(eL,:) + grad_comp(:,2)';
    elseif BE(i,4) == 2 || BE(i,4) == 3 || BE(i,4) == 4 % Wall edges
        uhat = un(eL,:);
        grad_comp = uhat'*BEnormal(i,:).*BElength(i);
        gradx(eL,:) = gradx(eL,:) + grad_comp(:,1)'; 
        grady(eL,:) = grady(eL,:) + grad_comp(:,2)';
    end
end

% Gradient divided by respective cell areas
gradx = gradx./area; grady = grady./area;

% Initializing Residual
R  = zeros(nelem, 4); smag = zeros(nelem,1);

% Interior edge flux
for i = 1:nIE
    n1 = IE(i,1); n2 = IE(i,2); eL = IE(i,3); eR = IE(i,4);

    % Edge midpoints
    emx = 1/2*(V(n1,1) + V(n2,1));    % Edge midpoint x-coordinates
    emy = 1/2*(V(n1,2) + V(n2,2));    % Edge midpoint y-coordinates
    em = [emx,emy];

    % x-xi
    dxL = em - cent(eL,:); dxR = em - cent(eR,:);

    % Face midpoint state
    gradiL = [gradx(eL,:)',grady(eL,:)']'; gradiR = [gradx(eR,:)',grady(eR,:)']';
    UL = un(eL,:) + (gradiL'*dxL')'; 
    UR = un(eR,:) + (gradiR'*dxR')';

    % Edge flux
    [Fi, smagi] = flux(UL,UR,IEnormal(i,:));

    % Edge residual
    R(eL,:) = R(eL,:) + Fi*IElength(i);
    R(eR,:) = R(eR,:) - Fi*IElength(i);

    % Propagation speed of disturbance
    smag(eL) = smag(eL) + abs(smagi)*IElength(i);    % Left cell 
    smag(eR) = smag(eR) + abs(smagi)*IElength(i);    % Right cell 
end

% Boundary edge flux
for i = 1:nBE
    n1 = BE(i,1); n2 = BE(i,2); eL = BE(i,3);

    % Edge midpoints
    emx = 1/2*(V(n1,1) + V(n2,1));    % Edge midpoint x-coordinates
    emy = 1/2*(V(n1,2) + V(n2,2));    % Edge midpoint y-coordinates
    em = [emx,emy];

    % x-xi
    dxL = em - cent(eL,:);

    % Face midpoint state
    gradiL = [gradx(eL,:)',grady(eL,:)']'; 
    UL = un(eL,:) + (gradiL'*dxL')'; 

    if BE(i,4) == 1 % Freestream edges
        [Fb, smagb] = flux(UL,uinf,BEnormal(i,:));          
    elseif BE(i,4) == 2 || BE(i,4) == 3 || BE(i,4) == 4     % Wall edges
        [Fb, smagb] = wallflux(UL,BEnormal(i,:));
    end

    % Edge residual
    R(eL,:) = R(eL,:) + Fb*BElength(i);

    % Propagation speed of disturbance
    smag(eL) = smag(eL) + abs(smagb)*BElength(i); % Boundary cell
end