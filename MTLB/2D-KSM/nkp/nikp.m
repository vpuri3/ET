function [X, Y] = nikp(A, B, mX, nX, mY, nY, tol, verbose)
% [X, Y] = nikp(A, B, mX, nX, mY, nY, tol, verbose)
%
% Given A and B, find X (mX-by-nX matrix) and Y (mY-by-nY matrix)
% to minimize || A (X kron Y) - B ||_F. Assume without loss of
% generality that ||X||_F = 1.
%
% Input:
%   A : m-by-(mX mY) matrix
%   B : m-by-(nX nY) matrix
%   mX, nX : integers > 0. Size of X.
%   mY, nY : integers > 0. Size of Y.
%   tol: scalar > 0 and << 1. Convergence tolerance.
%     Optional, default value is 1e-10.
%   verbose: 0 or 1. 1 means full verbosity.
%     Optional, default value is 0.
%
% Output:
%   X : mX-by-nX matrix. ||X||_F = 1.
%   Y : mY-by-nY matrix

if nargin < 6
  error('nikp error: at least 6 arguments required')
end

% Get sizes.
[mA, nA] = size(A);
[mB, nB] = size(B);

% Make sure dimensions match.
if mA ~= mB
  error('nikp error: size mismatch in A and B')
end
if nA ~= mX * mY
  error('nikp error: size mismatch in A and (X kron Y)')
end
if nX * nY ~= nB
  error('nikp error: size mismatch in (X kron Y) and B')
end

% Compensate for missing arguments.
if nargin < 7
  tol = 1e-10;
end
if nargin < 8
  verbose = 0;
end

% Convert from matrix problem to vector problem:
%   || A (X kron Y) - B ||_F = || A2 (x kron y) - b ||_2
% where
%   A2 = I_{nX^2} kron (I_{nY} kron A)(P_{mX, mX nY} kron I_{mY})
%   x = vec(X)
%   y = vec(Y)
%   b = vec(B)

% Compute w.
v = PerfShuf(mX, mX * nY);
lenw = mX * mY * nY;
w = zeros(1,lenw);
w(1 : mY : lenw) = (v-1) * mY + 1;
for j = 2 : mY
  w(j : mY : lenw) = w(1 : mY : lenw) + (j-1);
end

% Solve NIKPv problem.
[x, y] = nikpv('nikpmult', [mA, mX, nX, mY, nY, w], A, ...
    reshape(B, mB * nB, 1), mA*nX*nY, mX*nX, mY*nY, tol, verbose);
X = reshape(x, mX, nX);
Y = reshape(y, mY, nY);
