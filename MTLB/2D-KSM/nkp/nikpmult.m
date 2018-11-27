function y = nikpmult(mode, m, n, x, info, A)
% Required by nikp.m and LSQR.
%
% Compute products on
%   A2 = I_{nX} kron [ (I_{nY} kron A) I_{mX mY nY}(:,w) ]
%
% Input:
%   mode: 1 or 2. If mode = 1, compute y = A2 * x.
%     If mode = 2, compute y = A2' * x.
%   m, n : ignored
%   info : Required information packed into a vector:
%     info = [mA, mX, nX, mY, nY, w]
%     w is the permutation vector
%   A : mA-by-(mX mY) matrix

mA = info(1);
mX = info(2);
nX = info(3);
mY = info(4);
nY = info(5);
w = info(6 : mX*mY*nY+5);

if mode == 1
  X = reshape(x, mX*mY*nY, nX);
  Y = kronIA(nY, A, X(w,:));
  y = reshape(Y, mA*nX*nY, 1);
elseif mode == 2
  X = reshape(x, mA*nY, nX);
  Y(w,:) = kronIA(nY, A', X);
  y = reshape(Y, mX*mY*nY*nX, 1);
else
  error('nikpmult error: mode not recognized')
end
