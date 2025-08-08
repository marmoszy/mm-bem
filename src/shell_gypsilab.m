% Scattering from shell targets (air inside fluid)
% using GYPSILAB toolbox for Matlab
% MM 7.6.2025:  based on Gonzalez, Elavia 2020 

clear all; close all; clc; tic
run('../../gypsilab/addpathGypsilab.m')     % Gypsilab path

fname = '../msh/sphere-1.905-600.msh'; 
%fname = '../msh/YFT_swimbladder_origin.msh';
fname2 = '../msh/sphere-1.0-300.msh'; 
f0 = 38e3; 
th = 0;                % wave direction angle
c0 = 1480;             % water medium
rho0 = 1024; 
c1 = 1540;             % fluid medium (target)
rho1 = 1045;
c2 =  340;             % contaminated air medium
rho2 = 1.29; 

disp("Calculating ...");
k0 = 2*pi*f0/c0;          % sea water;
k1 = 2*pi*f0/c1;          % fish body
k2 = 2*pi*f0/c2;          % fish swimbladder
d = [cos(th*pi/180) sin(th*pi/180) 0];
rho01 = (rho0+rho1)/2;  
rho12 = (rho1+rho2)/2;

% mesh1
[vx,fx] = mshReadMsh(fname);
mesh1 = msh(vx,fx);        % mshSphere(600,0.01905);
N1 = size(mesh1.vtx,1);
sigma1 = dom(mesh1, 3);                         % quadrature
v1  = fem(mesh1, 'P1');                         % finite element
I1  = integral(sigma1, v1, v1);                 % identity
A11 = calderon2(sigma1, sigma1, v1, v1, k0, rho0); % scaled Calderons 
A11 = A11 - calderon2(sigma1, sigma1, v1, v1, k1, rho1);
A11 = A11 + [rho01*I1 zeros(size(I1));zeros(size(I1)) rho01*I1];

% mesh2
[vx2,fx2] = mshReadMsh(fname2);
mesh2 = msh(vx2,fx2);        % mshSphere(300, 0.01);
% mshWriteMsh('../msh/sphere-1.0-300.msh',mesh2);
N2 = size(mesh2.vtx,1);
sigma2 = dom(mesh2, 3);                         % quadrature
v2  = fem(mesh2, 'P1');                         % finite element
I2  = integral(sigma2, v2, v2);                 % identity
A22 = calderon2(sigma2, sigma2, v2, v2, k1, rho1);
A22 = A22 - calderon2(sigma2, sigma2, v2, v2, k2, rho2);
A22 = A22 + [rho12*I2 zeros(size(I2));zeros(size(I2)) rho12*I2];

% mesh1/mesh2 part 
A12 = -calderon2(sigma1, sigma2, v1, v2, k1, rho1);
A21 =  calderon2(sigma2, sigma1, v2, v1, k1, rho1);

% incident wave and its gradient
PW = @(X) exp(1i*k0*X*d');              
gradxPW{1} = @(X) 1i*k0*d(1).*PW(X);
gradxPW{2} = @(X) 1i*k0*d(2).*PW(X);
gradxPW{3} = @(X) 1i*k0*d(3).*PW(X);

% surface solution
f = integral(sigma1, v1, PW);               % incident wave traces
g = integral(sigma1, ntimes(v1), gradxPW); 
z = zeros(N2,1);
uu = [A11, A12; A21, A22]/rho0 \ [-f; g/rho0; z; z]; 

% far field solution
th = (0:359)' * pi/180;
r1 = [cos(th),sin(th),zeros(size(th))]; % far field cicle
[Sinf, Dinf] = potential2(r1, sigma1, v1, k0);
psc = rho0 * Dinf * uu(1:N1) + rho0.^2 * Sinf * uu(N1+1:2*N1); % Far field radiation
psc = psc/rho0;

% save, plot and print
s = [(0:359)' abs(psc)];
fid=fopen('../out/shell-gypsilab.txt','w');fprintf(fid,'%d\t%.6f\n',s');fclose(fid);
%!/usr/local/bin/gnuplot -c ../bin/polar.gp ../out/shell-gypsilab.txt
polarplot(th,max(-63,20*log10(abs(psc)))); rlim([-63 -20])
disp(['th0   = ' num2str(abs(psc(1))) ' (' num2str(20*log10(abs(psc(1)))) ')']);
disp(['th180 = ' num2str(abs(psc(length(psc)/2))) ' (' num2str(20*log10(abs(psc(length(psc)/2)))) ')']);
toc
