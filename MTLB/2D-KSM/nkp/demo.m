% Script for testing nikp.m

mX = 11;
nX = 23;
mY = 13;
nY = 17;
m = mX * nX + mY * nY + 10;

epsilon = 0;

A = 5 - 10 * rand(m, mX*mY);
X0 = 5 - 10 * rand(mX, nX);
Y0 = 5 - 10 * rand(mY, nY);
normx = norm(X0, 'fro');
X0 = X0 / normx;
Y0 = Y0 * normx;
B0 = 5 - 10 * rand(m, nX*nY);
B = A * kron(X0, Y0) + epsilon * B0;


[X, Y] = nikp(A, B, mX, nX, mY, nY);


fprintf('Residual = %e\n', norm(A * kron(X, Y) - B, 'fro'))
ex = min(norm(X - X0, 'fro'), norm(X + X0, 'fro'));
ey = min(norm(Y - Y0, 'fro'), norm(Y + Y0, 'fro'));
fprintf('Error X =  %e\n', ex)
fprintf('Error Y =  %e\n', ey)
