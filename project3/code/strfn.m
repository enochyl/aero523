function psi = strfn(mesh,un,order)
% PURPOSE: computes and plots the streamline of the given mesh, order of 
% method, and u-state
%
% INPUTS:
%   mesh  : input mesh structure, consisting of:
%     mesh.V  =  x,y coordinates of all nodes
%     mesh.E  =  triangular elements, as sequences of 3 nodes
%     mesh.IE =  interior faces and connectivities
%     mesh.BE =  boundary faces and connectivities
%   un    : input state
%   order : order of the method used
%
% OUTPUTS:
%   psi   : streamfunction
%

V = mesh.V;
E = mesh.E;
IE = mesh.IE;
BE = mesh.BE;

nnode = length( V);
nelem = length( E);
nIE   = length(IE);
nBE   = length(BE);

%% Interior edge normal and length
% Coordinates of interior nodes
ix1 = V(IE(1:nIE,1),1); ix2 = V(IE(1:nIE,2),1); 
iy1 = V(IE(1:nIE,1),2); iy2 = V(IE(1:nIE,2),2);

% Edge length
IElength = sqrt((ix1-ix2).^2 + (iy1-iy2).^2);

% Normal components
IEnormal(:,1) =  (iy2-iy1)./IElength;
IEnormal(:,2) = -(ix2-ix1)./IElength;

%% Boundary edge normal and length
% Coordinates of boundary nodes
bx1 = V(BE(1:nBE,1),1); bx2 = V(BE(1:nBE,2),1); 
by1 = V(BE(1:nBE,1),2); by2 = V(BE(1:nBE,2),2);

% Edge length
BElength = sqrt((bx1-bx2).^2 + (by1-by2).^2);

% Normal components
BEnormal(:,1) =  (by2-by1)./BElength;
BEnormal(:,2) = -(bx2-bx1)./BElength;

%% Cell area %%
% A = 1/2(x1y2 + x2y3 + x3y1 - x1y3 - x2y1 - x3y2)
temp  = 1/2*abs((V(E(:,1),1).*V(E(:,2),2) + V(E(:,2),1).*V(E(:,3),2)...
        + V(E(:,3),1).*V(E(:,1),2) - V(E(:,1),1).*V(E(:,3),2)...
        - V(E(:,2),1).*V(E(:,1),2) - V(E(:,3),1).*V(E(:,2),2)));
area = repmat(temp,1,4);

%% Centroids %%
cx = 1/3*(V(E(:,1),1)+V(E(:,2),1)+V(E(:,3),1)); % Centroid x-coordinates
cy = 1/3*(V(E(:,1),2)+V(E(:,2),2)+V(E(:,3),2)); % Centroid y-coordinates
cent = [cx,cy];
    
%% Streamfunction
psi = zeros(length(V),1);
[uinf] = init;  % Initialize u0
psicounter = zeros(length(V),1);    % Defining a vector that marks if a node's streamfunction value has been updated

for i = 1:nBE   % Setting the streamfunction marker of the boundary nodes on the main airfoil to be 1 such that their 0 values will not be changed later.
    if BE(i,4) == 2
        n1 = BE(i,1);
        n2 = BE(i,2);
        psicounter(n1) = 1;
        psicounter(n2) = 1;
    end
end

switch order
    case 1
        while sum(psicounter) < nnode   % Looping through all the nodes in the computational domain
            for i = 1:nIE   % Interior nodes
                eL = IE(i,3); eR = IE(i,4);
                uL = un(eL,:); uR = un(eR,:);
                [Fb, ~] = flux(uL,uR,IEnormal(i,:));
                
                if psicounter(IE(i,1)) == 1 && psicounter(IE(i,2)) == 0
                    psi(IE(i,2)) = psi(IE(i,1)) + Fb(1)*IElength(i);
                    psicounter(IE(i,2)) = psicounter(IE(i,2)) + 1;  % Updating the marker to prevent duplicated updates

                elseif psicounter(IE(i,1)) == 0 && psicounter(IE(i,2)) == 1
                    psi(IE(i,1)) = psi(IE(i,2)) - Fb(1)*IElength(i);
                    psicounter(IE(i,1)) = psicounter(IE(i,1)) + 1;  % Updating the marker to prevent duplicated updates
                end
            end
            for i = 1:nBE   % Boundary nodes
                eL = BE(i,3); uL = un(eL,:); 
                if BE(i,4) == 1
                    [Fb, ~] = flux(uL,uinf,BEnormal(i,:));              % Freestream edges
                elseif BE(i,4) == 2 || BE(i,4) == 3 || BE(i,4) == 4     % Wall edges
                    [Fb, ~, ~] = wallflux(uL,BEnormal(i,:));
                end
                
                if psicounter(BE(i,1)) == 0 && psicounter(BE(i,2)) == 1
                    psi(BE(i,1)) = psi(BE(i,2)) + Fb(1)*BElength(i);
                    psicounter(BE(i,1)) = psicounter(BE(i,1)) + 1;  % Updating the marker to prevent duplicated updates
                end

                if psicounter(BE(i,1)) == 1 && psicounter(BE(i,2)) == 0
                    psi(BE(i,2)) = psi(BE(i,1)) - Fb(1)*BElength(i);
                    psicounter(BE(i,2)) = psicounter(BE(i,2)) + 1;  % Updating the marker to prevent duplicated updates
                end
            end
        end
        
    case 2
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

        while sum(psicounter) < nnode
            for i = 1:nIE   % Interior nodes
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

                [Fb, ~] = flux(UL,UR,IEnormal(i,:));

                if psicounter(IE(i,1)) == 1 && psicounter(IE(i,2)) == 0
                    psi(IE(i,2)) = psi(IE(i,1)) + Fb(1)*IElength(i);
                    psicounter(IE(i,2)) = psicounter(IE(i,2)) + 1;  % Updating the marker to prevent duplicated updates

                elseif psicounter(IE(i,1)) == 0 && psicounter(IE(i,2)) == 1
                    psi(IE(i,1)) = psi(IE(i,2)) - Fb(1)*IElength(i);
                    psicounter(IE(i,1)) = psicounter(IE(i,1)) + 1;  % Updating the marker to prevent duplicated updates
                end
            end
            for i = 1:nBE   % Boundary nodes
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
                    [Fb, ~] = flux(UL,uinf,BEnormal(i,:));          
                elseif BE(i,4) == 2 || BE(i,4) == 3 || BE(i,4) == 4  % Wall edges
                    [Fb, ~] = wallflux(UL,BEnormal(i,:));
                end

                if psicounter(BE(i,1)) == 0 && psicounter(BE(i,2)) == 1
                    psi(BE(i,1)) = psi(BE(i,2)) + Fb(1)*BElength(i);
                    psicounter(BE(i,1)) = psicounter(BE(i,1)) + 1;  % Updating the marker to prevent duplicated updates
                end

                if psicounter(BE(i,1)) == 1 && psicounter(BE(i,2)) == 0
                    psi(BE(i,2)) = psi(BE(i,1)) - Fb(1)*BElength(i);
                    psicounter(BE(i,2)) = psicounter(BE(i,2)) + 1;  % Updating the marker to prevent duplicated updates
                end
            end    
        end
        
end


figure()
colormap(jet(2^5));
t = E';
t = [t; ones(1,size(t,2))];
pdeplot(V',[],t,'xydata',psi','xystyle','interp', 'colormap', 'jet','Contour','on','levels',12000); % last input for levels has to be changed depending on how many contours are desired.
axis equal; axis off;
axis([-0.3 0.9 -0.35 0.25]);
caxis([-0.1 0.1])








