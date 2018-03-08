function idsolve5(N)
% 2D Solver.
% Solve -(f*uxx + g*uyy) = pi^2*sin(pi*x)*sin(pi*y)
% domain [0,1]^2, homogeneous Dirichlet BC
% true solution u = sin(pi*x)*sin(pi*y)
if nargin == 0  N = 256; end
dx = 1/(N+1);
[X,Y] = ndgrid(0:dx:1);
T = sin(pi*X).*sin(pi*Y); T = T(2:end-1,2:end-1);
F = X.*ones(size(X)); F = F(2:end-1,2:end-1);
G = sin(exp(X.*Y)).*ones(size(X)); G = G(2:end-1,2:end-1);
R = T .* (F+G)*pi^2; % right hand side
r = reshape(R,N^2,1);

Ax = spdiags([-1*ones(N,1),2*ones(N,1),-1*ones(N,1)],-1:1,N,N);
Ax = Ax/(dx^2); % good ol' tridiagonal SPD

function [v] = afun(u)
    U = reshape(u,N,N);
    V = F.*(Ax*U) + G.*(U*Ax');
    v = reshape(V,N^2,1);

end      

Ax = full(Ax);
[Sx,Dx] = eig(Ax); DD = diag(Dx);
e = ones(N,1);
L = e*DD'+DD*e';
L = 1 ./L;

function [v] = pfun(u)
U = reshape(u,N,N);
V = Sx * (L .* (Sx' * U * Sx) ) * Sx';
v = reshape(V,N^2,1);
end

x = bicgstab(@afun,r,1e-10,500,@pfun);
X = reshape(x,N,N);
error = max(max(abs(X-T)))

end