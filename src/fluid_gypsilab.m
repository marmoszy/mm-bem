% Scattering from fluid targets 
% MM 6.6.2025 using GYPSILAB toolbox
% based on the equations from Gonzalez, Elavia 2020
clear all; close all; clc; tic
run('../../gypsilab/addpathGypsilab.m')  % Gypsilab path

fname = '../msh/sphere-1.905-600.msh'; 
%fname = '../msh/YFT_swimbladder_origin.msh';
f0 = 38e3;  th0 = 360;             % wave direction angle
c0 = 1480;  rho0 = 1024;           % water medium 
c1 = 1540;  rho1 = 1045;           % fluid medium (target)
oname = '../out/fluid-gypsilab.txt';

k0 = 2*pi*f0/c0;       % sea water;
k1 = 2*pi*f0/c1;       % fish body

disp("Assembling ...")
disp(fname) % input mesh
[vx,fx] = mshReadMsh(fname);
mesh = msh(vx,fx);       % mshSphere(600,0.01905); 
sigma = dom(mesh, 3);                    % quadrature
u  = fem(mesh, 'P1');                    % finite element

A0 = calderon2(sigma, sigma, u, u, k0, rho0);      % scaled Calderons 
A1 = calderon2(sigma, sigma, u, u, k1, rho1);
I  = integral(sigma,u,u);                          % identity
rho = (rho0+rho1)/2;
II = [rho*I zeros(size(I));zeros(size(I)) rho*I];  % scaled identity
[L,U] = lu(II + A0 - A1);

disp("Solving ..."); ss=[];
for th=0:th0
d = [cos(th*pi/180) sin(th*pi/180) 0];
% incident wave
PW = @(X) exp(1i*k0*X*d');              
gradxPW{1} = @(X) 1i*k0*d(1).*PW(X);
gradxPW{2} = @(X) 1i*k0*d(2).*PW(X);
gradxPW{3} = @(X) 1i*k0*d(3).*PW(X);

% surface solution
f = integral(sigma, u, PW);             % incident wave traces
g = integral(sigma, ntimes(u), gradxPW);  
uu  = U \ (L \ [-f; g/rho0]); %uu = (II + A0 - A1) \ [-f; g/rho0];
N = size(uu,1)/2;

% far field solution
th1 = (0:359)' * pi/180;
r1 = [cos(th1),sin(th1),zeros(size(th1))]; % far field cicle
xdoty   = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) exp(-1i*k0*xdoty(X,Y))/(4*pi);
gradxGinf{1} = @(X,Y) -1i*k0*X(:,1).*exp(-1i*k0*xdoty(X,Y))/(4*pi);
gradxGinf{2} = @(X,Y) -1i*k0*X(:,2).*exp(-1i*k0*xdoty(X,Y))/(4*pi);
gradxGinf{3} = @(X,Y) -1i*k0*X(:,3).*exp(-1i*k0*xdoty(X,Y)/(4*pi));
SL = integral(r1, sigma, Ginf, u);
DL = integral(r1, sigma, gradxGinf, ntimes(u)) ;
psc = rho0 * DL * uu(1:N) + rho0.^2 * SL * uu(N+1:2*N); % scattered

% save, plot and print
s = [(0:359)' abs(psc)]; mode=['w','a'];
fid=fopen(oname,mode((th~=0)+1));fprintf(fid,'%d\t%.6f\n',s');fprintf(fid,'\n\n');fclose(fid);
%!/usr/local/bin/gnuplot -c ../bin/polar.gp ../out/scat3-gypsilab.txt
polarplot(th1,max(-63,20*log10(abs(psc)))); rlim([-63 -20]);title(th);drawnow
i = mod(th+180,360)+1; % backscattering angle index
ss = [ss; s(i,:)];
disp([num2str(s(i,1)) ' ' num2str(s(i,2)) ]);
end; 
fid=fopen(strrep(oname,'.txt','-bsl.txt'),'w');fprintf(fid,'%d\t%.6f\n',ss');fclose(fid);
if size(ss,1)>1
    polarplot(ss(:,1)*pi/180,max(-63,20*log10(ss(:,2)))); rlim([-63 -20]);title('TS');
end
toc


