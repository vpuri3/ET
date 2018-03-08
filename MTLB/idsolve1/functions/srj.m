%% INPUTS
% Guess, True Solution Vectors, filename, # iterations.
%% OUTPUTS
% Solution vector after jacobi iteration, convergence plot.

function [A_out,srj_plot] = srj(A,R,T,V0,it)

load P6_350.srj
load P15_350.srj
load P6_150.srj
load P15_150.srj

Omega = P6_150;

M = length(Omega);
OneMinusOmega = 1 -  Omega;

n = size(A,1)-1;
h = 1/n;
hsq = h^2;
iteration = zeros(1,it);
residual = zeros(1,it);
res = max(max(abs(T-A)));
disp('STARTING SRJ ITERATION.')

%% k = iteration count + 1
k = 1 %Zeroth iteration
res = max(max(abs(T-A)));
residual(k) = res;
k = 2
%Some Constants
[X,Y] = ndgrid(0:h:1); %Domain is [0,1]x[0,1].
Rc = 3./(((X+1).^2).*(Y+1)) +3./(((Y+1).^2).*(X+1)); %some constants on RHS
laplaceV0 = V0
laplaceV0(2:n,2:n) = n^2*(V0(1:n-1,2:n)+V0(3:n+1,2:n)+V0(2:n,3:n+1)+V0(2:n,1:n-1)-4*V0(2:n,2:n));

while k <= it
    r = rem(k,M);
    if r == 0
        r = M;
    end
    B = 0.25*(A(1:n-1,2:n)+A(3:n+1,2:n)+A(2:n,3:n+1)+A(2:n,1:n-1)-hsq*R(2:n,2:n));
    A(2:n,2:n) = OneMinusOmega(r)*A(2:n,2:n) + Omega(r)*B;
    
    res = max(max(abs(T-A)));
    residual (k) = res;
    iteration (k) = k-1; % Iteration Count
    if rem(k,500) == 0
        iteration (k)
    end
    %% UNCOMMENT IF RHS HAS DEPENDENCE ON A
    R = 2*(V0.^3+3*V0.^2.*A-Rc)-laplaceV0; % update RHS
    k = k+1;
end

disp('Residual after SRJ Iteration is:')
res
A_out = A

%% PLOTTING
hold on
loglog(iteration,residual,'-b','LineWidth',1.2)
xlabel('Iteration Count','Fontsize',15)
ylabel('||r||_\infty','Fontsize',15)
legend('SRJ')
title('SRJ Iteration','Fontsize',15)
grid on
srj_plot = figure;

end