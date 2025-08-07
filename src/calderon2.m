function C = calderon2(sigma, sigma2, u, v, k, rho)
% Gonzalez scaled Calderon operator MM 6.6.2025

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
% grady Green kernel function --> G(x,y) = grady[exp(ik|x-y|)/|x-y|]
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Single layer --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy
S = 1/(4*pi) .* integral(sigma,sigma2,u,Gxy,v);
S = S + 1/(4*pi) .* regularize(sigma,sigma2,u,'[1/r]',v);
% Double layer --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy
D = 1/(4*pi) .* integral(sigma,sigma2,u,dyGxy,ntimes(v));
D = D + 1/(4*pi) .* regularize(sigma,sigma2,u,'grady[1/r]',ntimes(v));
% Double layer transpose --> \int_Sx \int_Sy psi(x)' dnx G(x,y) psi(y) dx dy
%Dt = D.';
Dt = 1/(4*pi) .* integral(sigma2,sigma,v,dyGxy,ntimes(u));
Dt = Dt + 1/(4*pi) .* regularize(sigma2,sigma,v,'grady[1/r]',ntimes(u));
Dt = Dt.';
% Hypersingular --> k^2 * \int_Sx \int_Sy n.psi(x) G(x,y) n.psi(y) dx dy
%                   - \int_Sx \int_Sy nxgrad(psi(x)) G(x,y) nxgrad(psi(y)) dx dy
H = 1/(4*pi) .* (k^2 * integral(sigma,sigma2,ntimes(u),Gxy,ntimes(v)) ...
    - integral(sigma,sigma2,nxgrad(u),Gxy,nxgrad(v)));
H = H + 1/(4*pi) .* (k^2 * regularize(sigma,sigma2,ntimes(u),'[1/r]',ntimes(v)) ...
    - regularize(sigma,sigma2,nxgrad(u),'[1/r]',nxgrad(v)));
%C = [-D S ; -H Dt];               % Calderon operator
C = [rho*D rho*rho*S ; -H -rho*Dt];  % scaled2 Calderon operator
