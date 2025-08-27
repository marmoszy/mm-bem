% scattering from soft target
% this script requires the GYPSILAB toolbox for Matlab!           |
% MM 3.8.2025
clear all; close all; clc; tic
run('../../gypsilab/addpathGypsilab.m')  % Gypsilab path

% params
%fname = '../msh/sphere-1.905-600.msh'; 
fname = '../msh/YFT_swimbladder_origin.msh';
Th0=-90:10:90; Ph0=0:10:360;       % wave direction angles
f0 = 38e3; 
c0 = 1480;  
%oname='../out/soft-gypsilab.txt';
oname = '../out/YFT_swimbladder_origin-gypsilab.txt';

disp("Assembling ...");
k0 = 2*pi*f0/c0;

% input mesh
[vx,fx] = mshReadMsh(fname);
mesh = msh(vx,fx);       % mshSphere(600,0.01905); 
sigma = dom(mesh,3);     % Domain gss=3
u = fem(mesh,'P0');      % Finite elements typ='P0'|'P1';
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k0);
S = 1/(4*pi) .* integral(sigma,sigma,u,Gxy,u,1e-5);   % tol=1e-5
%Sr  = 1/(4*pi) .* regularize(sigma,sigma,u,'[1/r]',u);
%S = S + Sr;
[L,U] = lu(S);

disp("Solving ..."); ss=[];
for ph=Ph0

% preparation for far field solution 
th1 = pi/180 .* (1:360)';
phr=ph*pi/180;
rotph=[cos(phr) 0 sin(phr); 0 1 0; -sin(phr) 0 cos(phr)];
r1 = (rotph*[cos(th1),sin(th1),zeros(size(th1))]')';
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) 1/(4*pi) .* exp(-1i*k0*xdoty(X,Y));
SL = integral(r1,sigma,Ginf,u,1e-5);

for th=Th0
% Incident wave
    d = (rotph*[cos(th*pi/180) sin(th*pi/180) 0]')';
    PW = @(X) exp(1i*k0*X*d'); 

% Surface solution 
    f = integral(sigma,u,PW);
    lambda  = U \ (L \ -f);  % lambda = S \ -f;

% far field solution
    psc = SL * lambda;  

% save, plot and print
    s = [(0:359)' abs(psc)]; mode=['w','a'];
    fid=fopen(oname,mode((th~=0)+1));fprintf(fid,'%d\t%.6f\n',s');fprintf(fid,'\n');fclose(fid);
    %!/usr/local/bin/gnuplot -c ../bin/polar.gp ../out/soft-gypsilab.txt
    polarplot(th1,max(-63,20*log10(abs(psc)))); rlim([-63 -20]);
    title([num2str(th),' ',num2str(ph)]);drawnow
    i = mod(th+180,360)+1; % backscattering angle index
    ss = [ss; s(i,1) ph s(i,2)];
    %disp([num2str(s(i,1)) ' ' num2str(s(i,2)) ]);
end;
end;

fid2=fopen(strrep(oname,'.txt','-bsl.txt'),'w');
fprintf(fid2,sprintf('# soft_gypsilab # %dx%d\n',length(Th0),length(Ph0)));
fprintf(fid2,'%d\t%d\t%.6f\n',ss');
fclose(fid2);
if size(ss,1)>1
    polarplot(ss(:,1)*pi/180,max(-63,20*log10(ss(:,3)))); rlim([-63 -20]);title('TS');
end
toc

