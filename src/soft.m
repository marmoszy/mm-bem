% plane wave scattering from soft sphere in salt water
% input: msh ascii 2.2 file
% output: 1-deg scattering pattern in txt file
% MM 26.7.2025            |
clear all; close all; clc; tic;

% input parameters
fname = '../msh/sphere-1.905-600.msh'; 
f = 38e3; c = 1480; k = 2*pi*f/c; X0  = [1 0 0];

% reading vtx,elt from mesh file
fid = fopen(fname,'r');
str = fgets(fid); str = fgets(fid); str = fgets(fid); str = fgets(fid);
Nvtx = str2double(fgets(fid));
vtx  = zeros(Nvtx,3);
for i = 1:Nvtx,tmp = str2num(fgets(fid)); vtx(i,:) = tmp(2:4); end
str = fgets(fid); str = fgets(fid);
Nelt = str2double(fgets(fid));
elt  = zeros(Nelt,3);
for i = 1:Nelt,tmp = str2num(fgets(fid)); elt(i,1:3) = tmp(6:8);end
fclose(fid);

% space 
X = (vtx(elt(:,1),:)+vtx(elt(:,2),:)+vtx(elt(:,3),:))/3;
%E1 = vtx(elt(:,2),:) - vtx(elt(:,1),:);
%E2 = vtx(elt(:,3),:) - vtx(elt(:,1),:);
%Mu = diag(0.5*sqrt(sum(cross(E1,E2,2).^2,2))); 
Mu = 1 ;

% Helmholtz single layer operator
Nx = size(X,1);
S = zeros(Nx,Nx);
for j = 1:Nx
    r = sqrt((X(:,1)-X(j,1)).^2+(X(:,2)-X(j,2)).^2+(X(:,3)-X(j,3)).^2);
    S(:,j) =  1/(4*pi) * exp(1i*k*r)./r;
    S(r<1e-12,j) = 0 + 1i*k/(4*pi);
end

% surface solution
LHS =  S * Mu;
RHS = - exp(1i*k*X*X0');
lambda = LHS \ RHS; %[Lh,Uh] = lu(LHS); lambda  = Uh \ (Lh \ RHS);

% far field solution
th = pi/180 .* (1:360)';
r1 = [cos(th),sin(th),zeros(size(th))];
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
G = zeros(size(r1,1),Nx);
for j = 1:Nx, G(:,j) = 1/(4*pi) .* exp(-1i*k*xdoty(r1,X(j,:))); end
SL = G * Mu;  
psc = SL * lambda;    

% save result
s = [(0:359)' abs(psc)];
fid=fopen('../txt/soft-m.txt','w');fprintf(fid,'%d\t%.6f\n',s');fclose(fid);
!/usr/local/bin/gnuplot -p -c ../txt/polar.gp ../txt/soft-m.txt
disp(['th0   = ' num2str(abs(psc(1))) ' (' num2str(20*log10(abs(psc(1)))) ')']);
disp(['th180 = ' num2str(abs(psc(length(psc)/2))) ' (' num2str(20*log10(abs(psc(length(psc)/2)))) ')']);
toc

