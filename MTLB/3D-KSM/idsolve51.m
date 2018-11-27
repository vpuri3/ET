function idsolve51(N)
% 3D Solver.
% Solve -(uxx + uyy+ uzz) = 3*pi^2*sin(pi*x)*sin(pi*y)
% domain [0,1]^2, homogeneous Dirichlet BC
% true solution u = sin(pi*x)*sin(pi*y)*sin(pi*z)
if nargin == 0  N = 32; end
dx = 1/(N+1);
[X,Y,Z] = ndgrid(0:dx:1);
T = sin(pi*X).*sin(pi*Y).*sin(pi*Z); % true solution
T = T(2:end-1,2:end-1,2:end-1);
R = 3*pi^2*T; % right hand side
r = reshape(R,N^3,1);

function [v] = afun(u)
    % produces the affect of Laplace operator on vector u.
    U = reshape(u,N,N,N);
    UU = zeros(N+2,N+2,N+2);
    UU(2:N+1,2:N+1,2:N+1) = U;
    UU(2:N+1,2:N+1,2:N+1) = (-1/dx^2)*(UU(1:N,2:N+1,2:N+1)+UU(3:N+2,2:N+1,2:N+1) + UU(2:N+1,1:N,2:N+1)+UU(2:N+1,3:N+2,2:N+1) + UU(2:N+1,2:N+1,1:N)+UU(2:N+1,2:N+1,3:N+2) - 6*U);
    V = UU(2:N+1,2:N+1,2:N+1);
    v = reshape(V,N^3,1);
end

x = cgs(@afun,r);
X = reshape(x,N,N,N);
error = max(max(max( abs(X-T) )))

end