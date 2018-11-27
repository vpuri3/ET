% solve cxx*u_xx + cyy*u_yy + c1*u = r

n = 100; %Number of grid points per dimension.
h = 1/n; %Grid spacing
it_jacobi = 1000;

[X,Y] = ndgrid(0:h:1);
%Domain is [0,1]x[0,1]. %First index is for x (top to bottom).
%Second index is for y axis (left to right).
hsq = h^2;
nsq = n^2;

cxx = ones(size(X));
cyy = ones(size(X));%Y.^3 + 1;
c1 = exp(X+Y);
T = exp(X.*Y);%exp(X.*Y);
Txx = Y.*Y;
Tyy = X.*X;
R = cxx.*Txx + cyy.*Tyy + c1.*T;

A = zeros(n+1);
A(:,1) = T(:,1);
A(:,n+1) = T(:,n+1);
A(1,:) = T(1,:);
A(n+1,:) = T(n+1,:);
% RHS

res = max(max(abs(T-A)));
iteration = zeros(1,it_jacobi+1);
residual = zeros(1,it_jacobi+1);
residual (1) = res;

onebyfour = 0.25;

disp('STARTING JACOBI ITERATION.')
omega = 1;

count = 1;
while count <= it_jacobi
    S = R(2:n,2:n) - nsq*(cxx(2:n,2:n).*(A(1:n-1,2:n)+A(3:n+1,2:n))+cyy(2:n,2:n).*(A(2:n,1:n-1)+A(2:n,3:n+1)));
    S = S ./ (c1(2:n,2:n)-2*nsq*(cxx(2:n,2:n)+cyy(2:n,2:n)));
    A(2:n,2:n) = (1-omega)*A(2:n,2:n) + omega*S;
    res = max(max(abs(T-A)));%norm(T-A,2);
    residual (count+1) = res;
    iteration (count+1) = count;
    count = count+1;
end

hold on
plot(iteration(1:it_jacobi+1),residual(1:it_jacobi+1),'-','LineWidth',2);

legend('SOR')
xlabel('Iteration Count','Fontsize',12)
ylabel('|| r ||','Fontsize',12)
grid on
