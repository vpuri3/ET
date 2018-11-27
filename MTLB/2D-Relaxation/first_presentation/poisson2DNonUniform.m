clc
%Input Parameters
nx = 585;
ny = 280;
%om = [1] %Array of relaxation factors.

Omega = om;
it = 1e4; %Number of iterations.
hx = 1/nx; %Grid spacing
hy = 1/ny;

%N_ef = pi/(asin(sqrt(2*(sin(pi/(2*n)))^2)));

%Setting stuff up.
syms x y;
pointsX = 0:hx:1;
pointsY = 0:hy:1;
[X,Y] = ndgrid(pointsX,pointsY);
T = zeros([nx+1,ny+1]);
R = zeros([nx+1,ny+1]);

%True Solution
true = -exp(x*y);
T = -exp(X.*Y);
t = norm(T,inf);

%First index is for x (up to down). Second index is for y axis (left to right).
%Filling in Dirchlet boundary conditions. Domain is [0,1]x[1,0].
A = zeros([nx+1,ny+1]);
A(:,1) = T(:,1);
A(:,ny+1) = T(:,ny+1);
A(1,:) = T(1,:);
A(nx+1,:) = T(nx+1,:);

%Source Term = RHS
rhs = true*(x^2+y^2);
R = T.*(X.^2+Y.^2);


iteration = zeros(1,it);
residual = ones(1,it);
res = norm(T-A, inf);

%Computing constants outside loop.
M = length(om);
%hsq = h^2;

OneMinusOmega = 1 -  Omega;
%iteration count
k = 1
while k < it
    r = rem(k,M);
    if r == 0
        r = M;
    end
    B = 0.25*(A(1:n-1,2:n)+A(3:n+1,2:n)+A(2:n,3:n+1)+A(2:n,1:n-1)-hsq*R(2:n,2:n));
    A(2:n,2:n) = OneMinusOmega(r)*A(2:n,2:n) + Omega(r)*B;
    if rem(k,1e1) == 0
    res = norm((T-A),inf)/t;
    end
    residual (k) = res;
    R = A;
    iteration (k) = k;
    k = k + 1
end

res

%plotting
close all
loglog(iteration,residual,'r','LineWidth',1)
xlabel('Iteration Count','Fontsize',12)
ylabel('Relative Error','Fontsize',12)
legend('Convergence P')
title('256x256')
grid on