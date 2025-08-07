% Scattering from fluid targets
% using GYPSILAB toolbox for Matlab 
% MM 6.6.2025 - based on the equations from Gonzalez, Elavia 2020
clear all; close all; clc; tic
run('../../gypsilab/addpathGypsilab.m')           % Gypsilab path

fname = '../msh/sphere-1.905-600.msh'; 
%fname = '../msh/YFT_swimbladder_origin.msh';
f0 = 38e3; 
th = 0;                % wave direction angle
c0 = 1480;             % water medium
rho0 = 1024; 
c1 = 1540;             % fluid medium (target)
rho1 = 1045;

k0 = 2*pi*f0/c0;       % sea water;
k1 = 2*pi*f0/c1;       % fish body
d = [cos(th*pi/180) sin(th*pi/180) 0];

% input mesh
[vx,fx] = mshReadMsh(fname);
mesh = msh(vx,fx);       % mshSphere(600,0.01905); 
sigma = dom(mesh, 3);                    % quadrature
u  = fem(mesh, 'P1');                    % finite element

A0 = calderon2(sigma, sigma, u, u, k0, rho0);      % scaled Calderons 
A1 = calderon2(sigma, sigma, u, u, k1, rho1);
I  = integral(sigma,u,u);                          % identity
rho = (rho0+rho1)/2;
II = [rho*I zeros(size(I));zeros(size(I)) rho*I];  % scaled identity

% incident wave
PW = @(X) exp(1i*k0*X*d');              
gradxPW{1} = @(X) 1i*k0*d(1).*PW(X);
gradxPW{2} = @(X) 1i*k0*d(2).*PW(X);
gradxPW{3} = @(X) 1i*k0*d(3).*PW(X);

% surface solution
f = integral(sigma, u, PW);             % incident wave traces
g = integral(sigma, ntimes(u), gradxPW);  
uu = (II + A0 - A1) \ [-f; g/rho0];
N = size(uu,1)/2;

% far field solution
th = (0:359)' * pi/180;
r1 = [cos(th),sin(th),zeros(size(th))]; % far field cicle
xdoty   = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) exp(-1i*k0*xdoty(X,Y))/(4*pi);
gradxGinf{1} = @(X,Y) -1i*k0*X(:,1).*exp(-1i*k0*xdoty(X,Y))/(4*pi);
gradxGinf{2} = @(X,Y) -1i*k0*X(:,2).*exp(-1i*k0*xdoty(X,Y))/(4*pi);
gradxGinf{3} = @(X,Y) -1i*k0*X(:,3).*exp(-1i*k0*xdoty(X,Y)/(4*pi));
Sinf = integral(r1, sigma, Ginf, u);
Dinf = integral(r1, sigma, gradxGinf, ntimes(u)) ;
psc = rho0 * Dinf * uu(1:N) + rho0.^2 * Sinf * uu(N+1:2*N); % Far field radiation

% plot, print and save 
s = [(0:359)' abs(psc)];
fid=fopen('../out/fluid-gypsilab.txt','w');fprintf(fid,'%d\t%.6f\n',s');fclose(fid);
%!/usr/local/bin/gnuplot -p -c ../bin/polar.gp ../out/fluid-gypsilab.txt
polarplot(th,max(-63,20*log10(abs(psc)))); rlim([-63 -20])
disp(['th0   = ' num2str(abs(psc(1))) ' (' num2str(20*log10(abs(psc(1)))) ')']);
disp(['th180 = ' num2str(abs(psc(length(psc)/2))) ' (' num2str(20*log10(abs(psc(length(psc)/2)))) ')']);
toc;



