% notes and stuff

function [v] = afun(u)
% produces the affect of Laplace operator on vector u.
U = reshape(u,N,N,N);
V = zeros(N,N,N);
W = zeros(N,N,N);
for k = 1:N
    V(:,:,k) = Ax*U(:,:,k) + U(:,:,k)*Ax';
end
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