%% This is 2D Tensor Product Solver with scaled coefficients.
% 3uxx + 4uyy = -r
clear;
n = 8;
h = 1/(n+1);
[X,Y] = ndgrid(linspace(h,1-h,n)); % same as linspace(0,1,n)
T = sin(pi*X).*sin(pi*Y);
R = (pi^2)*T*7;

I = speye(n);
A = spdiags([-1*ones(n,1),2*ones(n,1),-1*ones(n,1)],-1:1,n,n);
A = A / (h^2);
Ax = 3*full(A); Ay = 4*full(A);

[Vx,Dx] = eig(Ax);
[Vy,Dy] = eig(Ay);
DDx = diag(Dx); DDy = diag(Dy);

% A = Iy x Ax + Ay x Ix = (Vy x Vx)(Ix x Dx + Dy x Ix) (inv(Vy) x inv(Vx))
% (P x Q) u = Q*U*P'
e = ones(n,1);
L = e*DDy'+DDx*e';
L = 1 ./L;
U = Vx * (L .* (Vx' * R * Vy) ) * Vy';
error = max(max(abs(T-U)))