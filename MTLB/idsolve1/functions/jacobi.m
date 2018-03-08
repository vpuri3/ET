function [A_out,jacobi_plot] = jacobi(A,R,T,it)
%% INPUTS
% Guess Solution Vector, True Solution Vector, # iterations.
%% OUTPUTS
% Solution vector after jacobi iteration, convergence plot.

n = size(A,1)-1;
h = 1/n;
hsq = h^2;
iteration = zeros(1,it);
residual = zeros(1,it);
disp('STARTING JACOBI ITERATION.')

[X,Y] = ndgrid(0:h:1);
Rc = 3./(((X+1).^2).*(Y+1)) +3./(((Y+1).^2).*(X+1)); %some constants on RHS
%% k = iteration count + 1
k = 1 %Zeroth iteration
res = max(max(abs(T-A)));
residual(k) = res;
k = 2

while k <= it
    A(2:n,2:n) = 0.25*(A(1:n-1,2:n)+A(3:n+1,2:n)+A(2:n,3:n+1)+A(2:n,1:n-1)-hsq*R(2:n,2:n));
    res = max(max(abs(T-A)));
    residual (k) = res;
    iteration (k) = k-1; % Iteration Count
    if rem(k,300) == 0
        iteration (k)
    end
    %% UNCOMMENT IF RHS HAS DEPENDENCE ON A
    R = 2*(A.^3-Rc); % update RHS
    k = k+1;
end

disp('Residual after Jacobi Iteration is:')
res
A_out = A;

%% PLOTTING
hold on;
loglog(iteration,residual,'-k','LineWidth',1.2);
xlabel('Iteration Count','Fontsize',15);
ylabel('||r||_\infty','Fontsize',15);
legend('Jacobi');
title('Jacobi Iteration','Fontsize',15);
grid on;
jacobi_plot = figure;

end