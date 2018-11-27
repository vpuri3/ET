function [x, y] = nikpv(Amult, iw, rw, b, m, n1, n2, tol, verbose)
% [x, y] = nikpv(Amult, iw, rw, b, m, n1, n2, tol, verbose)
%
% Solve the NIKPv problem
%   || A (x kron y) - b ||_2
% where A is m-by-(n1*n2), x is an n1-vector, y is an n2-vector,
% b is an m-vector. Amult is the name of a function which computes
% products on A. iw and rw are passed into Amult (for example,
% sparse storage of A).

% 1. Find w to minimize || A w - b ||_2
% 2. Find x and y to minimize || x kron y - w ||_2

% Solve w = A \ b.
tol2 = tol / 100;
w = lsqr(m, n1*n2, Amult, iw, rw, b, 0, tol2, tol2, 1e8, ...
   n1*n2, verbose);

% Solve min || w - x kron y ||.
B = reshape(w, n2, n1)';

[U, S, V] = svd(B);
x = U(:,1);
y = V(:,1) * S(1,1);
