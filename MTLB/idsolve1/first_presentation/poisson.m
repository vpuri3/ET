clc
clearvars -except om15_150 om15_200 om6_150

n = 256; %Number of grid points per dimension.
om = 1; %Array of relaxation factors.
Omega = om;

it = 10e3; %Number of iterations.
h = 1/n; %Grid spacing

N_ef = pi/(2*(asin(sqrt(2*(sin(pi/(2*n)))^2))));

points = 0:h:1;
[X,Y] = ndgrid(points);
T = zeros(n+1);
R = zeros(n+1);

%True Solution
T = 1./(X+1)+1./(Y+1);
t = max(max(T)); %infinity norm of True Solution

%First index is for x (top to bottom). Second index is for y axis (left to right).
%Filling in Dirchlet boundary conditions. Domain is [0,1]x[1,0].
A = zeros(n+1);
A(:,1) = T(:,1);
A(:,n+1) = T(:,n+1);
A(1,:) = T(1,:);
A(n+1,:) = T(n+1,:);

%Source Term = RHS
Rc = 3./(((X+1).^2).*(Y+1)) +3./(((Y+1).^2).*(X+1)); %some constants in right hand side
R = 2*(T.^3-Rc);

iteration = zeros(1,it);
residual = zeros(1,it);
res = max(max(abs(T-A)));

%Computing constants outside loop.
M = length(Omega);
hsq = h^2;
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
    
    if rem(k,1) == 0
    res = max(max(abs(T-A)));
    end
    residual (k) = res;
    %R = 2*(A.^3-Rc);
    iteration (k) = k;
    k = k + 1
end

res

%plotting
%close all
loglog(iteration,residual,'k','LineWidth',1)
xlabel('Iteration Count','Fontsize',15)
ylabel('||r||_(inf)','Fontsize',15)
%legend('SRJ P15, N150','SRJ P15, N200','Jacobi')
title('Dirichlet, Source phi^3+..., True 1/(x+1)+1/(y+1), 256x256','Fontsize',15)
grid on
hold on