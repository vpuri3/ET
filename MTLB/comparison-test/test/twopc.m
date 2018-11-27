load P6_100.srj
Omega = P6_100;
n = 100;
dmsz = 5;
h = 2*dmsz/n;
it_jacobi = 200;
it_srj = length(Omega);
it_linear = 3;
points = (-1*dmsz):h:dmsz;
N_ef = 0.5*pi/(asin(sqrt((2/3)*(sin(pi/(2*n)))^2))); % Dirichlet Boundary Conditions;

disp(['Domain: [' num2str(points(1)) ',' num2str(points(end)) ']^3' ])
disp(['Grid Spacing: ' num2str(h) ', ' num2str(n+1) ' points per dimension'])
disp(['Jacobi Iterations: ' num2str(it_jacobi)])
disp(['SRJ Iterations: ' num2str(it_srj)])
disp(['Newton Raphson Steps: ' num2str(it_linear)])
disp(['Effective number of SRJ points: ' num2str(N_ef)])
%First index is for x (top to bottom). Second index is for y axis (left to right).
%Third index is for z axis.

[X,Y,Z] = ndgrid(points);
% Getting Bowen York Curvature.
% 1,2,3 --> x,y,z respectively
c = 1.168642873; % position of black holes
p = [0, 0.3331917498, 0];
M_plus = 0.453;
M_minus = 0.453;

c_plus = [c,0,0];             % position of black hole
c_minus = -1*c_plus;
r_plus = sqrt(Y.^2+Z.^2+(X-c_plus(1)).^2);
r_minus = sqrt(Y.^2+Z.^2+(X-c_minus(1)).^2);

p_plus = zeros(n+1,n+1,n+1,3); 
p_plus(:,:,:,1) = ones(n+1,n+1,n+1) *p(1) ;
p_plus(:,:,:,2) = ones(n+1,n+1,n+1) *p(2) ;
p_plus(:,:,:,3) = ones(n+1,n+1,n+1) *p(3) ;

p_minus = -1*p_plus;

n_plus = zeros(n+1,n+1,n+1,3);
n_plus(:,:,:,1) = (X-c_plus(1)) ./ r_plus;
n_plus(:,:,:,2) = Y ./ r_plus;
n_plus(:,:,:,3) = Z ./ r_plus;

n_minus = zeros(n+1,n+1,n+1,3);
n_minus(:,:,:,1) = (X-c_minus(1)) ./ r_minus;
n_minus(:,:,:,2) = Y ./ r_minus;
n_minus(:,:,:,3) = Z ./ r_minus;

dot_plus = zeros(n+1,n+1,n+1);
dot_minus = zeros(n+1,n+1,n+1);
for w = 1:3
    dot_plus  = dot_plus  + n_plus (:,:,:,w).*p_plus (:,:,:,w);
    dot_minus = dot_minus + n_minus(:,:,:,w).*p_minus(:,:,:,w);
end

Asq = zeros(n+1,n+1,n+1);
Aij = zeros(n+1,n+1,n+1);
for m = 1:3
for s = 1:3
    Aij = + 1.5 * (p_plus(:,:,:,m) .* n_plus(:,:,:,s) + p_plus(:,:,:,s) .* n_plus(:,:,:,m) + dot_plus .* n_plus(:,:,:,m) .* n_plus(:,:,:,s)) ./ (r_plus.^2);
    Aij = Aij + 1.5*(p_minus(:,:,:,m) .* n_minus(:,:,:,s) + p_minus(:,:,:,s) .* n_minus(:,:,:,m) + dot_minus .* n_minus(:,:,:,m) .* n_minus(:,:,:,s)) ./ (r_minus.^2);
    if (m == s)
        Aij = Aij - +1.5 * (dot_plus./(r_plus.^2) + dot_minus./(r_minus.^2));
    end
    Asq = Asq + Aij .* Aij;
end
end

C = 1+0.5*(M_plus./r_plus)+0.5*(M_minus./r_minus);
% gxx is psi^4 and psi = 1 + par_m_plus/rplus + par_m_minus/rminus + u
h5disp('gxx.xyz.h5','/ADMBASE::gxx it=0 tl=0 rl=0')
T = h5read('gxx.xyz.h5','/ADMBASE::gxx it=0 tl=0 rl=0');
T = T.^0.25;
T = T-C;

disp('Test Parameters')
disp(['M_plus= ' num2str(M_plus) ', M_minus= ' num2str(M_minus)])
disp([ 'Black Hole Plus Location: (' num2str(c_plus(1)) ',0,0)' ])
disp([ 'Black Hole Minus Location: (' num2str(c_minus(1)) ',0,0)' ])
disp([ 'Black Hole Plus Momentum: (' num2str(p(1)) ',' num2str(p(2)) ',' num2str(p(3)) ')' ])