% plane wave scattering from soft sphere in salt water
% input: msh ascii 2.2 file
% output: 1-deg scattering pattern in txt file
% MM 26.7.2025            |
clear all; close all; clc; tic;

% input parameters
fname = '../msh/sphere-1.905-600.msh'; 
%fname = '../msh/YFT_swimbladder_origin.msh';
th= 270;       % plane wave direction angle in xy plane
f = 38e3;    % its frequency
c = 1480;    % sound speed in water

k = 2*pi*f/c; 
d = [cos(th*pi/180) sin(th*pi/180) 0];

% reading vx,fx from mesh file
fid = fopen(fname,'r');
fgets(fid); fgets(fid); fgets(fid); fgets(fid); % skip 4 lines
N = str2double(fgets(fid));
vx  = zeros(N,3);
for i = 1:N,tmp = str2num(fgets(fid)); vx(i,:) = tmp(2:4); end
fgets(fid); fgets(fid); % skip 2 $lines
M = str2double(fgets(fid));
fx  = zeros(M,3);
for i = 1:M,tmp = str2num(fgets(fid)); fx(i,1:3) = tmp(6:8);end
fclose(fid);

% P0 space - triangle centers
X = (vx(fx(:,1),:)+vx(fx(:,2),:)+vx(fx(:,3),:))/3;
% P1 space - three points from a triangular element
%X = zeros(3*M,3); w = [4 1 1; 1 4 1; 1 1 4]/6;
%for j = 1:3
%    X((j:3:3*M)',:)=w(j,1)*vx(fx(:,1),:)+w(j,2)*vx(fx(:,2),:)+w(j,3)*vx(fx(:,3),:);
%end

% Helmholtz single layer operator
Nx = size(X,1);
S = zeros(Nx,Nx);
for j = 1:Nx
    r = sqrt((X(:,1)-X(j,1)).^2+(X(:,2)-X(j,2)).^2+(X(:,3)-X(j,3)).^2);
    S(:,j) =  exp(1i*k*r)./r /(4*pi);
    S(r<1e-12,j) = 0 + 1i*k/(4*pi);
end

%incident field
f = exp(1i*k*X*d');

% surface solution
u = S \ -f; 

% far field solution
th = pi/180 .* (1:360)';
r1 = [cos(th),sin(th),zeros(size(th))];
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
SL = zeros(size(r1,1),Nx);
for j = 1:Nx, SL(:,j) = exp(-1i*k*xdoty(r1,X(j,:)))/(4*pi) ; end
psc = SL * u;    

% plot, print and save 
s = [(0:359)' abs(psc)];
fid=fopen('../out/soft-m.txt','w');fprintf(fid,'%d\t%.6f\n',s');fclose(fid);
%!/usr/local/bin/gnuplot -p -c ../bin/polar.gp ../out/soft-m.txt
polarplot(th,20*log10(abs(psc))); rlim([-63 -20]);
disp(['th0   = ' num2str(abs(psc(1))) ' (' num2str(20*log10(abs(psc(1)))) ')']);
disp(['th180 = ' num2str(abs(psc(length(psc)/2))) ' (' num2str(20*log10(abs(psc(length(psc)/2)))) ')']);
toc
