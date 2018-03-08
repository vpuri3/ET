clc
%Input Parameters
n = 50; %Number of grid points per dimension.

om = [1]; %Array of relaxation factors.

Omega = om;
it = 0.1e4; %Number of iterations.
h = 1/n; %Grid size

%N_ef = pi/(asin(sqrt(2*(sin(pi/(2*n)))^2)));

%Setting stuff up.
syms x y z;
points = 0:h:1;
[X,Y,Z] = ndgrid(points,points,points);
T = zeros(n+1);
R = zeros(n+1);

%True Solution
true = x+y;
T = X+Y;
t = sqrt(sum(sum(sum(T.*T))));

%First index is for x (up to down). Second index is for y axis (left to right).
%Filling in Dirchlet boundary conditions. Domain is [1,0,0]x[0,1,0]x[0,0,1].
A = zeros([n+1,n+1,n+1]);
A(:,:,1) = T(:,:,1);
A(:,:,n+1) = T(:,:,n+1);
A(:,1,:) = T(:,1,:);
A(:,n+1,:) = T(:,n+1,:);
A(1,:,:) = T(1,:,:);
A(n+1,:,:) = T(n+1,:,:);

%Source Term = RHS
rhs = 0;
R = 0*ones([n+1,n+1,n+1]);

iteration = zeros(1,it);
residual = ones(1,it);
res = sqrt(sum(sum(sum((T-A).*(T-A)))));

%Computing constants outside loop.

M = length(om);
hsq = h^2;
OneMinusOmega = 1 -  Omega;
%iteration count
k = 1
fac = 1/4;
while k < it
    r = rem(k,M);
    if r == 0
        r = M;
    end
    B = fac*(A(1:n-1,2:n,2:n)+A(3:n+1,2:n,2:n)+A(2:n,1:n-1,2:n)+A(2:n,3:n+1,2:n)+A(2:n,2:n,1:n-1)+A(2:n,2:n,3:n+1)-hsq*R(2:n,2:n,2:n));
    A(2:n,2:n,2:n) = OneMinusOmega(r)*A(2:n,2:n,2:n) + Omega(r)*B;
    if rem(k,3e2) == 0
    res = sqrt(sum(sum(sum((T-A).*(T-A)))))/t;
    end
    residual (k) = res;
    R = A;
    iteration (k) = k;
    k = k + 1
end

res

%plotting
loglog(iteration,residual)
xlabel('Iteration Count','Fontsize',12)
ylabel('Relative Error','Fontsize',12)
grid on
axis square