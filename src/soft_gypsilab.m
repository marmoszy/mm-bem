% scattering from soft target
% this script requires the GYPSILAB toolbox for Matlab!           |
% MM 3.8.2025
clear all; close all; clc; tic
run('../../gypsilab/addpathGypsilab.m')  % Gypsilab path

% params
fname = '../msh/sphere-1.905-600.msh'; 
%fname = '../msh/YFT_swimbladder_origin.msh';
th=0; 
f = 38e3; 
c = 1480;  

k = 2*pi*f/c;  %d  = [0 -1 0];
d = [cos(th*pi/180) sin(th*pi/180) 0];

% input mesh
[vx,fx] = mshReadMsh(fname);
mesh = msh(vx,fx);       % mshSphere(600,0.01905); 
sigma = dom(mesh,3);     % Domain gss=3
u = fem(mesh,'P0');      % Finite elements typ='P0'|'P1';

% Incident wave
PW = @(X) exp(1i*k*X*d'); 

% Surface solution 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
S = 1/(4*pi) .* integral(sigma,sigma,u,Gxy,u,1e-5);   % tol=1e-5
Sr  = 1/(4*pi) .* regularize(sigma,sigma,u,'[1/r]',u);
S = S + Sr;
f = integral(sigma,u,PW);
%[Lh,Uh] = lu(S);lambda  = Uh \ (Lh \ f); 
lambda = S \ -f;

% far field solution
th = pi/180 .* (1:360)';
r1 = [cos(th),sin(th),zeros(size(th))];
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) 1/(4*pi) .* exp(-1i*k*xdoty(X,Y));
SL = integral(r1,sigma,Ginf,u);
psc = SL * lambda;  

% plot, print and save 
s = [(0:359)' abs(psc)];
fid=fopen('../out/soft-gypsilab.txt','w');fprintf(fid,'%d\t%.6f\n',s');fclose(fid);
%!/usr/local/bin/gnuplot -p -c ../bin/polar.gp ../out/soft-gypsilab.txt
polarplot(th,20*log10(abs(psc))); rlim([-63 -20])
disp(['th0   = ' num2str(abs(psc(1))) ' (' num2str(20*log10(abs(psc(1)))) ')']);
disp(['th180 = ' num2str(abs(psc(length(psc)/2))) ' (' num2str(20*log10(abs(psc(length(psc)/2)))) ')']);
toc;

