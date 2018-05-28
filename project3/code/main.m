clear
close all
clc

% PURPOSE: main driver

%% User input
mesh_size = input('Enter mesh size (0/1/2/3): ');
order = input('Order of simulation (1/2): ');
inputcheck(mesh_size,order);
mesh = mesh(mesh_size);

%% Mesh
plotmesh(mesh);

%% Flowfield state solver
tic
switch order
    case 1
        [un,iter,Rmax,cl,cd,cl_conv,cd_conv,xpos,cpd] = solvefv1(mesh);
    case 2
        [un,iter,Rmax,cl,cd,cl_conv,cd_conv,xpos,cpd] = solvefv2(mesh);
    otherwise
        error('unsupported');
end
toc
%% Mach contour
plotsoln(mesh, un, 'mach', order);

%% Pressure contour
plotsoln(mesh, un, 'pressure',order);

%% Rmax convergence
figure()
n = 1:1:iter;
plot(n,Rmax)
xlabel('Iteration'); ylabel('log_{10}R_{max}');
grid on

%% Cl convergence
figure()
plot(n,cl_conv)
xlabel('Iteration'); ylabel('c_l');
grid on
cl_tol = cl_conv(end);  % total cl from the converged solution

%% Cd convergence
figure()
plot(n,cd_conv)
xlabel('Iteration'); ylabel('c_d');
grid on
cd_tol = cd_conv(end);  % total cd from the converged solution

%% Surface pressure coefficient
figure()
plot(xpos,cpd,'*')
xlabel('x position'); ylabel('c_p');
set(gca,'Ydir','reverse');
grid on

%% Streamline
psi = strfn(mesh,un,order);

%% Mass flowrate
[m_slat,m_flap] = flowrate(mesh,psi);




