function test1d(N)
% this is to test conjugate gradient with variable coefficients
% Solve -f*uxx = pi^2*sin(pi*x)
% domain [0,1]^1, homogeneous Dirichlet BC
% true solution u = sin(pi*x)
if nargin == 0
    N = 32;
end
dx = 1/(N+1);
x = 0:dx:1;
t = sin(pi*x); t = t(2:N+1)';
f = x.*ones(size(x)); f = f(2:N+1)';
r = pi^2 * t;

r = diag(f) * r;

Ax = spdiags([-1*ones(N,1),2*ones(N,1),-1*ones(N,1)],-1:1,N,N);
Ax = Ax/(dx^2); % SPD tridiagonal SPD

Ax = diag(f) * Ax;

s = bicgstab(Ax,r,1e-10,100);
error = max( abs(s-t) )

end