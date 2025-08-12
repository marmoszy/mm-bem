% Scattering from shell targets (i.e. air inside fluid)
% MM 7.6.2025 using gypsilab toolbox, calderon2.m, potential2.m
% equations based on Gonzalez, Elavia 2020 

clear all; close all; clc; tic;
run('../../gypsilab/addpathGypsilab.m')     % Gypsilab path

% input parameters
fname1 = '../msh/sphere-1.905-600.msh'; 
fname2 = '../msh/sphere-1.0-300.msh'; 
f0 = 38e3; th0 = 360;     % wave frequency and direction angle
c0 = 1480; rho0 = 1024;   % water medium 
c1 = 1540; rho1 = 1045;   % fluid medium (target)
c2 =  340; rho2 = 1.29;   % contaminated air medium
oname = '../out/shell-gypsilab.txt';

disp("Assembling ...");
k0 = 2*pi*f0/c0;          % sea water;
k1 = 2*pi*f0/c1;          % fish body
k2 = 2*pi*f0/c2;          % fish swimbladder
rho01 = (rho0+rho1)/2;  
rho12 = (rho1+rho2)/2;

disp(fname1);             % mesh1
[vx,fx] = mshReadMsh(fname1);
mesh1 = msh(vx,fx);       % mshSphere(600,0.01905);
% mshWriteMsh('mesh1.msh',mesh1);
N1 = size(mesh1.vtx,1);
sigma1 = dom(mesh1, 3);   % quadrature
v1  = fem(mesh1, 'P1');   % finite element
I1  = integral(sigma1, v1, v1); % identity
A11 = calderon2(sigma1, sigma1, v1, v1, k0, rho0); % scaled Calderons 
A11 = A11 - calderon2(sigma1, sigma1, v1, v1, k1, rho1);
A11 = A11 + [rho01*I1 zeros(size(I1));zeros(size(I1)) rho01*I1];

disp(fname2);              % mesh2
[vx2,fx2] = mshReadMsh(fname2);
mesh2 = msh(vx2,fx2);     % mshSphere(300, 0.01);
% mshWriteMsh('mesh2.msh',mesh2);
N2 = size(mesh2.vtx,1);
sigma2 = dom(mesh2, 3);     % quadrature
v2  = fem(mesh2, 'P1');     % finite element
I2  = integral(sigma2, v2, v2);   % identity
A22 = calderon2(sigma2, sigma2, v2, v2, k1, rho1);
A22 = A22 - calderon2(sigma2, sigma2, v2, v2, k2, rho2);
A22 = A22 + [rho12*I2 zeros(size(I2));zeros(size(I2)) rho12*I2];

% mesh1/mesh2 part 
A12 = -calderon2(sigma1, sigma2, v1, v2, k1, rho1);
A21 =  calderon2(sigma2, sigma1, v2, v1, k1, rho1);

disp("Solving ..."); ss=[];
for th=0:th0
% incident wave and its gradient
d = [cos(th*pi/180) sin(th*pi/180) 0];
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
th1 = (0:359)' * pi/180;
r1 = [cos(th1),sin(th1),zeros(size(th1))];  % far field circle
[SL, DL] = potential2(r1, sigma1, v1, k0);
psc = rho0 * DL * uu(1:N1) + rho0.^2 * SL * uu(N1+1:2*N1); % scattered
psc = psc/rho0;

% save, plot and print
s = [(0:359)' abs(psc)]; mode=['w','a'];
fid=fopen(oname,mode((th~=0)+1));fprintf(fid,'%d\t%.6f\n',s');fprintf(fid,'\n');fclose(fid);
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