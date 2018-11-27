%% This is 3D code with Symmetric Domain.

clear
close
n = 64; %Number of grid points per dimension.
dmsz = 0.65;
h = 2*dmsz/n; %Grid spacing
it_jacobi = 10;
it_srj = 1000;
it_linear = 3; % # of linearization steps.
points = (-1*dmsz):h:dmsz;
N_ef = 0.5*pi/(asin(sqrt((2/3)*(sin(pi/(2*n)))^2))) % Dirichlet Boundary Conditions;

disp(['Domain: [' num2str(points(1)) ',' num2str(points(end)) ']^3' ])
disp(['Grid Spacing: ' num2str(h) ', ' num2str(n) ' points per dimension'])
disp(['Jacobi Iterations: ' num2str(it_jacobi)])
disp(['SRJ Iterations: ' num2str(it_srj)])
disp(['Linearization Steps: ' num2str(it_linear)])
disp(['Effective number of SRJ points: ' num2str(N_ef)])

%% Initializing Stuff
[X,Y,Z] = ndgrid(points);
%First index is for x (top to bottom).
%Second index is for y (left to right).
%Third index is for z.
hsq = h^2;
% True Solution
T = 1./(X+1)+1./(Y+1);
%Dirchlet boundary conditions.
A = zeros(n+1,n+1,n+1);
A(1,:,:) = T(1,:,:);
A(n+1,:,:) = T(n+1,:,:);
A(:,1,:) = T(:,1,:);
A(:,n+1,:) = T(:,n+1,:);
A(:,:,1) = T(:,:,1);
A(:,:,n+1) = T(:,:,n+1);

% RHS
Rc = 3./(((X+1).^2).*(Y+1)) +3./(((Y+1).^2).*(X+1));
R = 2*(A.^3-Rc);

res = max(max(max(abs(T-A))));
iteration = zeros(1,it_jacobi+it_linear*it_srj+1);
residual = zeros(1,it_jacobi+it_linear*it_srj+1);
residual (1) = res;

onebysix = 1/6;

%% Jacobi Iteration on nonlinear PDE.
disp('Running Jacobi.')

count = 1;
while count <= it_jacobi
    A(2:n,2:n,2:n) = onebysix*(A(1:n-1,2:n,2:n)+A(3:n+1,2:n,2:n)+A(2:n,1:n-1,2:n)+A(2:n,3:n+1,2:n)+A(2:n,2:n,1:n-1)+A(2:n,2:n,3:n+1)-hsq*R(2:n,2:n,2:n));
    res = max(max(max(abs(T-A))));
    residual (count+1) = res;
    iteration (count+1) = count;
    %% UNCOMMENT IF RHS HAS DEPENDENCE ON A
    R = 2*(A.^3-Rc); % update RHS
    count = count+1;
end


%% Plotting
semilogy(iteration(1:it_jacobi+1),residual(1:it_jacobi+1),'--k','LineWidth',1.2,'DisplayName','Jacobi')
hold on

V0 = A;


%% Linearization and SRJ
load P6_64.srj
load P6_100.srj
load P6_300.srj

Omega = P6_64;
M = length(Omega);
OneMinusOmega = 1 -  Omega;

count=1; % Count
while count <= it_linear
    disp(['Newton Raphson ' num2str(count)])
    deltaV = zeros(n+1,n+1,n+1); % guess
    laplaceV0 = V0;
    laplaceV0(2:n,2:n,2:n) = (1/hsq)*(V0(1:n-1,2:n,2:n)+V0(3:n+1,2:n,2:n)+V0(2:n,1:n-1,2:n)+V0(2:n,3:n+1,2:n)+V0(2:n,2:n,1:n-1)+V0(2:n,2:n,3:n+1)-6*V0(2:n,2:n,2:n));
    Rv = 2*(V0.^3+3*V0.^2.*deltaV-Rc)-laplaceV0;
    Tv = T - V0;
    
    m = it_jacobi+it_srj*(count-1)+1;
    k=1;
    while k <= it_srj
        r = rem(k,M);
        if r == 0
            r = M;
        end
        B = onebysix*(deltaV(1:n-1,2:n,2:n)+deltaV(3:n+1,2:n,2:n)+deltaV(2:n,1:n-1,2:n)+deltaV(2:n,3:n+1,2:n)+deltaV(2:n,2:n,1:n-1)+deltaV(2:n,2:n,3:n+1)-hsq*Rv(2:n,2:n,2:n));
        deltaV(2:n,2:n,2:n) = OneMinusOmega(r)*deltaV(2:n,2:n,2:n) + Omega(r)*B;

        Diff = abs(Tv(2:n,2:n,2:n)-deltaV(2:n,2:n,2:n));
        res = max(max(max(Diff)));
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
