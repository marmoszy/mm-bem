function [Sinf,Dinf] = potentials(r, sigma, v, k)
% Green kernel function and its gradient
xdoty   = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) exp(-1i*k*xdoty(X,Y))/(4*pi);
gradxGinf{1} = @(X,Y) -1i*k*X(:,1).*exp(-1i*k*xdoty(X,Y))/(4*pi);
gradxGinf{2} = @(X,Y) -1i*k*X(:,2).*exp(-1i*k*xdoty(X,Y))/(4*pi);
gradxGinf{3} = @(X,Y) -1i*k*X(:,3).*exp(-1i*k*xdoty(X,Y)/(4*pi));

% far field solution
Sinf = integral(r, sigma, Ginf, v);
Dinf = integral(r, sigma, gradxGinf, ntimes(v)) ;