%% This is 2D code.

clear

n = 256; %Number of grid points per dimension.
h = 1/n; %Grid spacing
it_jacobi = 100;
it_srj = 1000;
it_linear = 7; % # of linearization steps.

N_ef = pi/(2*(asin(sqrt(2*(sin(pi/(2*n)))^2))));

%% Initializing Stuff
[X,Y] = ndgrid(0:h:1);
%Domain is [0,1]x[0,1]. %First index is for x (top to bottom).
%Second index is for y axis (left to right).
hsq = h^2;
nsq = n^2;
% True Solution
T = 1./(X+1)+1./(Y+1);
%infinity norm of True Solution.
%Dirchlet boundary conditions.
A = zeros(n+1);
A(:,1) = T(:,1);
A(:,n+1) = T(:,n+1);
A(1,:) = T(1,:);
A(n+1,:) = T(n+1,:);
% RHS
Rc = 3./(((X+1).^2).*(Y+1)) +3./(((Y+1).^2).*(X+1));
R = 2*(A.^3-Rc);

res = max(max(abs(T-A)));
iteration = zeros(1,it_jacobi+it_linear*it_srj+1);
residual = zeros(1,it_jacobi+it_linear*it_srj+1);
residual (1) = res;

onebyfour = 0.25;

%% Jacobi Iteration on nonlinear PDE.
disp('STARTING JACOBI ITERATION.')

count = 1;
while count <= it_jacobi
    A(2:n,2:n) = 0.25*(A(1:n-1,2:n)+A(3:n+1,2:n)+A(2:n,1:n-1)+A(2:n,3:n+1)-hsq*R(2:n,2:n));
    res = max(max(abs(T-A)));
    residual (count+1) = res;
    iteration (count+1) = count;
    %% UNCOMMENT IF RHS HAS DEPENDENCE ON A
    R = 2*(A.^3-Rc); % update RHS
    count = count+1;
end


%% Plotting
hold on
semilogy(iteration(1:it_jacobi+1),residual(1:it_jacobi+1),'--k','LineWidth',1.2,'DisplayName','Jacobi')


V0 = A;

%% Linearization and SRJ
load P6_150.srj
load P15_150.srj

Omega = P15_150;
M = length(Omega);
OneMinusOmega = 1 -  Omega;

disp('Starting with Linearization')

count=1; % Count
while count <= it_linear
    disp(['Linearization Step ' num2str(count)])
    deltaV = zeros(n+1); % guess
    laplaceV0 = V0;
    laplaceV0(2:n,2:n) = nsq*(V0(1:n-1,2:n)+V0(3:n+1,2:n)+V0(2:n,3:n+1)+V0(2:n,1:n-1)-4*V0(2:n,2:n));
    Rv = 2*(V0.^3+3*V0.^2.*deltaV-Rc)-laplaceV0;
    Tv = T - V0;
    
    m = it_jacobi+it_srj*(count-1)+1;
    k=1;
    while k <= it_srj
        r = rem(k,M);
        if r == 0
            r = M;
        end
        B = 0.25*(deltaV(1:n-1,2:n)+deltaV(3:n+1,2:n)+deltaV(2:n,3:n+1)+deltaV(2:n,1:n-1)-hsq*Rv(2:n,2:n));
        deltaV(2:n,2:n) = OneMinusOmega(r)*deltaV(2:n,2:n) + Omega(r)*B;

        Diff = abs(Tv(2:n,2:n)-deltaV(2:n,2:n));
        res = max(max(Diff));
        residual (m+k) = res;
        iteration (m+k) = m+k; % Iteration Count
        %% UNCOMMENT IF RHS HAS DEPENDENCE ON A
        Rv = 2*(V0.^3+3*V0.^2.*deltaV-Rc)-laplaceV0; % update RHS
        k = k+1;
    end
    
    %Plotting
    semilogy(iteration(m:m+k-1),residual(m:m+k-1),'-','Linewidth',1.2,'DisplayName',['SRJ ' num2str(count)])
    
    V0 = V0 + deltaV;
    count = count+1;
end

legend('show')
xlabel('Iteration Count','Fontsize',15)
ylabel('|| r ||_\infty','Fontsize',15)
title('Linearized SRJ','Fontsize',15)
grid on
