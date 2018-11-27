%% This is 2D Poisson Tensor Product Solver. uxx + uyy = -r
clear; % clc;
n = 1024;
h = 1/(n+1);
[X,Y] = ndgrid(linspace(h,1-h,n)); % same as linspace(0,1,n)
R = sin(pi*X).*sin(pi*Y);
T = (1/(2*pi*pi)).*R;

Ix = speye(n);
%Ix = spdiags([1*ones(n,1)],0,n,n);
Ax = spdiags([-1*ones(n,1),2*ones(n,1),-1*ones(n,1)],-1:1,n,n);
Ax = Ax / (h*h);
Ax = full(Ax);
Ay = Ax; Iy = Ix; % Rx = Ix(2:end-1,:); Ry = Rx;

[V,D] = eig(Ax);
DD = diag(D);
num = 1:n;
for c = num; % Orthagonalize eigenvectors for repeated eigenvalues.
    s = find(DD(c+1:end) == D(c));
    if ~isempty(s)
        s = s+c; s = [c,s];
        P = V(:,s);
        P = orth(P);
        V(:,s) = P;
        num(s) = [];
    end
end

% A = Iy x Ax + Ay x Ix = (Vy x Vx)(Ix x Dx + Dy x Ix) (inv(Vy) x inv(Vx))
e = ones(n,1);
L = e*DD'+DD*e';
L = 1 ./L;
U = V * (L .* (V' * R * V) ) * V';

error = max(max(abs(T-U)))