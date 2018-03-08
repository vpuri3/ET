function idsolve52(N)
% 3D Solver with scaled coefficients.
% Solve -(f*uxx + g*uyy+ h*uzz) = pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)(f+g+h)
% domain [0,1]^2, homogeneous Dirichlet BC
% true solution u = sin(pi*x)*sin(pi*y)*sin(pi*z)

if nargin == 0  N = 256; end
dx = 1/(N+1);
[X,Y,Z] = ndgrid(0:dx:1);
T = sin(pi*X).*sin(pi*Y).*sin(pi*Z); T = T(2:end-1,2:end-1,2:end-1);
F = ones(size(X)); F = F(2:end-1,2:end-1,2:end-1);
G = ones(size(X)); G = G(2:end-1,2:end-1,2:end-1);
H = ones(size(X)); H = H(2:end-1,2:end-1,2:end-1);
R = (pi^2)* T .* (F+G+H);
r = reshape(R,N^3,1);

UU = zeros(N+2,N+2,N+2); % embedd U into UU. Add BC to UU.

function [v] = afun(u)
    % produces the affect of scaled Laplace operator on vector u.
    U = reshape(u,N,N,N);
    UU(2:N+1,2:N+1,2:N+1) = U;
    V = (-1/dx^2)* F.* (UU(1:N,2:N+1,2:N+1)+UU(3:N+2,2:N+1,2:N+1) - 2*U);
    V = V + (-1/dx^2)* G.* (UU(2:N+1,1:N,2:N+1)+UU(2:N+1,3:N+2,2:N+1) - 2*U);
    V = V + (-1/dx^2)* H.* (UU(2:N+1,2:N+1,1:N)+UU(2:N+1,2:N+1,3:N+2) - 2*U);
    v = reshape(V,N^3,1);
end

function [v] = pfun(u)
   U = reshape(u,N,N,N); 
   FF = (F+G+H)/3;
   
end

s = bicgstab(@afun,r,1e-10,500);
S = reshape(s,N,N,N);
error = max(max(max( abs(S-T) )))

end