function test1d(N)
if nargin == 0
    N = 8;
end
dx = 1/(N+1);
x = 0:dx:1;
t = sin(pi*x); t = t(2:N+1)';
p = ones(size(x)); p = p(2:N+1)';
r = pi^2 * t;

%r = diag(f) * r;

D = spdiags([-ones(N,1),1*ones(N,1)],-1:0,N,N);
%Ax = Ax/(dx); % SPD tridiagonal SPD

A = D' * diag(p) * D;
A = full(A)


end