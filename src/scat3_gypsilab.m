% Scattering from three spheres (two inside a larger one) 
% MM 22.7.2025 using gypsilab toolbox,calderon2.m, potential2.m
% self developed equations!

clear all; close all; clc; tic
run('../../gypsilab/addpathGypsilab.m') % Gypsilab path

% input parameters
%a1 = 0.01905; N1 = 600;                 % sphere 1 params
%a2 = 0.005;   N2 = 200;                 % sphere 2 params
%a3 = 0.005;   N3 = 200;                 % sphere 3 params
fname1 = '../msh/sphere-1.905-600.msh'; 
fname2 = '../msh/sphere-0.5-200p.msh';
fname3 = '../msh/sphere-0.5-200n.msh'; 
f0 = 38e3; th0 = 360;    % wave freqency and direction
c0 = 1480; rho0 = 1024;  % sea water
c1 = 1540; rho1 = 1045;  % fish body
c2 =  340; rho2 = 1.29;  % fish swimbladder
c3 =  340; rho3 = 1.29;  % fish swimbladder
oname  = '../out/scat3-gypsilab.txt';

disp("Assembling ...");
k0 = 2*pi*f0/c0; k1 = 2*pi*f0/c1; k2 = 2*pi*f0/c2; k3 = 2*pi*f0/c3;
rho01 = (rho0+rho1)/2; rho12 = (rho1+rho2)/2; rho13 = (rho1+rho3)/2;

%mesh1  = mshSphere(N1, a1); 
%mshWriteMsh('mesh1.msh',mesh1);
disp(fname1);
[vx1,fx1] = mshReadMsh(fname1);
mesh1 = msh(vx1,fx1);           
N1 = size(mesh1.vtx,1);
sigma1 = dom(mesh1, 3);                         % quadrature
v1  = fem(mesh1, 'P1');                         % finite element
I1  = integral(sigma1, v1, v1);                 % identity
A11 = calderon2(sigma1, sigma1, v1, v1, k0, rho0); % scaled Calderon 
A11 = A11 - calderon2(sigma1, sigma1, v1, v1, k1, rho1);
A11 = A11 + [rho01*I1 zeros(size(I1));zeros(size(I1)) rho01*I1];

%mesh2 = mshSphere(N2, a2);
%mesh2 = mesh2.translate([0.01,0,0]);            % mesh sphere
%mshWriteMsh('mesh2.msh',mesh2);
disp(fname2);
[vx2,fx2] = mshReadMsh(fname2);
mesh2 = msh(vx2,fx2);           
N2 = size(mesh2.vtx,1);
sigma2 = dom(mesh2, 3);                         % quadrature
v2  = fem(mesh2, 'P1');                         % finite element
I2  = integral(sigma2, v2, v2);                 % identity
A22 = calderon2(sigma2, sigma2, v2, v2, k1, rho1);
A22 = A22 - calderon2(sigma2, sigma2, v2, v2, k2, rho2);
A22 = A22 + [rho12*I2 zeros(size(I2));zeros(size(I2)) rho12*I2];
A12 = -calderon2(sigma1, sigma2, v1, v2, k1, rho1);
A21 =  calderon2(sigma2, sigma1, v2, v1, k1, rho1);

%mesh3  = mshSphere(N3, a3);
%mesh3 = mesh3.translate([-0.01,0,0]);           % mesh sphere
%mshWriteMsh('mesh3.msh',mesh3);
disp(fname3);
[vx3,fx3] = mshReadMsh(fname3);
mesh3 = msh(vx3,fx3);            
N3 = size(mesh3.vtx,1);
sigma3 = dom(mesh3, 3);                         % quadrature
v3  = fem(mesh3, 'P1');                         % finite element
I3  = integral(sigma3, v3, v3);                 % identity
A33 = calderon2(sigma3, sigma3, v3, v3, k1, rho1);
A33 = A33 - calderon2(sigma3, sigma3, v3, v3, k3, rho3);
A33 = A33 + [rho13*I3 zeros(size(I3));zeros(size(I3)) rho13*I3];
A13 = -calderon2(sigma1, sigma3, v1, v3, k1, rho1);
A31 =  calderon2(sigma3, sigma1, v3, v1, k1, rho1);
Z23 = 0*calderon2(sigma2, sigma3, v2, v3, k1, rho1);
Z32 = 0*calderon2(sigma3, sigma2, v3, v2, k1, rho1);

disp("Solving ..."); ss=[];
for th=0:th0
% incident wave
d = [cos(th*pi/180) sin(th*pi/180) 0];
PW = @(X) exp(1i*k0*X*d');              
gradxPW{1} = @(X) 1i*k0*d(1).*PW(X);
gradxPW{2} = @(X) 1i*k0*d(2).*PW(X);
gradxPW{3} = @(X) 1i*k0*d(3).*PW(X);

% surface solution
f = integral(sigma1, v1, PW);               % incident traces
g = integral(sigma1, ntimes(v1), gradxPW); 
z2 = zeros(N2,1); z3 = zeros(N3,1);
uu = [A11, A12, A13; A21, A22, Z23; A31, Z32, A33] \ [-f; g/rho0; z2; z2; z3; z3];

% far field solution
th1 = (0:359)' * pi/180;
r1 = [cos(th1),sin(th1),zeros(size(th1))];
[SL, DL] = potential2(r1, sigma1, v1, k0);
psc = rho0 * DL * uu(1:N1) + rho0.^2 * SL * uu(N1+1:2*N1); % scattered field

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