%% This is 2D Tensor Product Solver with variable coefficients.
% f(x,y)uxx + g(x,y)uyy = -r;
n = 8;
h = 1/(n+1);
[X,Y] = ndgrid(linspace(h,1-h,n));
T = sin(pi*X).*sin(pi*Y);
F = 3*ones(size(X)); G = 4*ones(size(X));
%F = X; G = zeros(size(X));
R = (pi^2)*(F+G).*T;

I = speye(n);
A = spdiags([-1*ones(n,1),2*ones(n,1),-1*ones(n,1)],-1:1,n,n);
A = A / (h^2);
Ax = full(A); Ay = Ax;

[Vx,Dx] = eig(Ax);
[Vy,Dy] = eig(Ay);
DDx = diag(Dx); DDy = diag(Dy);
% A2D =  Df(Iy x Ax) + Dg(Ay x Ix)
% A2D = (Vy x Vx)(Ix x Dx + Dy x Ix) (inv(Vy) x inv(Vx))
% (P x Q) u = Q*U*P'
e = ones(n,1);
L = F.*(DDx*e') + G.*(e*DDy');
[U,S,V] = svd(L);
%L = 1 ./L;
%U = Vx * (L .* (Vx' * R * Vy) ) * Vy';
%U = (Vx' * R * Vy);
%error = max(max(abs(T-U)))