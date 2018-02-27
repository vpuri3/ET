%% Section 1
n = 50;
dmsz = 5;
h = 2*dmsz/n;
it_linear = 5;
it_jacobi = 50;
points = (-1*dmsz):h:dmsz;

disp(['Domain: [' num2str(points(1)) ',' num2str(points(end)) ']^3' ])
disp(['Grid Spacing: ' num2str(h) ', ' num2str(n+1) ' total points per dimension'])
disp(['Jacobi Iterations: ' num2str(it_jacobi)])
disp(['Linearization Steps: ' num2str(it_linear)])
%First index is for x (top to bottom). Second index is for y axis (left to right).
%Third index is for z axis.

[X,Y,Z] = ndgrid(points);
%% suffix 1 --> plus, suffix 2 --> minus
%syms c p r1 r2
p = 0.3331917498;
c = 1.168642873;
c1 = [c,0,0];
c2 = -1*c1;
p1 = [0,p,0];
p2 = -1*p1;
s1 = [0,0,0];
s2 = [0,0,0];
M1 = 0.453;
M2 = 0.453;
r1 = sqrt(Y.^2+Z.^2+(X-c1(1)).^2);
r2 = sqrt(Y.^2+Z.^2+(X-c2(1)).^2);
C = 1+(M1./r1)+(M2./r2);

% gxx is psi^4 and psi = 1 + par_m_plus/rplus + par_m_minus/rminus + u
h5disp('gxx.xyz.h5','/ADMBASE::gxx it=0 tl=0 rl=0')
T = h5read('gxx.xyz.h5','/ADMBASE::gxx it=0 tl=0 rl=0');
T = T.^0.25;
T = T-C;

% Components of Ap1
A11 = (3*p*Y.*(((c - X).*(c - X))./(r1.^2) - 1))./(2*r1.^3) ;
A12 = -(3*((p*(c - X))./r1 + (p*Y.*Y.*(c - X))./(r1.^3)))./(2*r1.^2) ;
A13 = -(3*p*Z.*Y.*(c - X))./(2*r1.^5) ;
A21 = -(3*((p*(c - X))./r1 + (p*Y.^2.*(c - X))./(r1.^3)))./(2*r1.^2) ;
A22 = (3*((Y*p)./r1 + (p*Y)./r1 + (p*Y.*((Y.^2)./(r1.^2) - 1))./r1))./(2*r1.^2) ;
A23 = (3*((Z*p)./r1 + (p*Z.*Y.^2)./(r1.^3)))./(2*r1.^2) ;
A31 = -(3*p*Y.*Z.*(c - X))./(2*r1.^5) ;
A32 = (3*((p*Z)./r1 + (p*Y.*Y.*Z)./(r1.^3)))./(2*r1.^2) ;
A33 = (3*p*Y.*((Z.*Z)./(r1.^2) - 1))./(2*r1.^3) ;

% Adding components of Ap2
A11 = A11 + -(3*p*Y.*(((c + X).*(c + X))./(r2.*r2) - 1))./(2*r2.^3) ;
A12 = A12 + -(3*((p*(C + X))./r2 + (p*Y.*Y.*(c + X))./(r2.^3)))./(2*r2.^2) ;
A13 = A13 + -(3*p*Z.*Y.*(c + X))./(2*r2.^5) ;
A21 = A21 + -(3*((p*(c + X))./r2 + (p*Y.^2.*(c + X))./(r2.^3)))./(2*r2.^2) ;
A22 = A22 + -(3*((Y*p)./r2 + (p*Y)./r2 + (p*Y.*((Y.^2)./(r2.^2) - 1))./r2))./(2*r2.^2) ;
A23 = A23 + -(3*((Z*p)./r2 + (p*Z.*Y.^2)./(r2.^3)))./(2*r2.^2) ;
A31 = A31 + -(3*p*Y.*Z.*(c + X))./(2*r2.^5) ;
A32 = A32 + -(3*((p*Z)./r2 + (p*Y.*Y.*Z)./(r2.^3)))./(2*r2.^2) ;
A33 = A33 + -(3*p*Y.*((Z.*Z)./(r2.^2) - 1))./(2*r2.^3) ;

Asq = A11.^2 + A12.^2 + A13.^2 + A21.^2 + A22.^2 + A23.^2 + A31.^2 + A32.^2 + A33.^2; 

%% Stuff
%v1 = cross(s1,n1);
%vv1 = n1'*v1;
%As1 = (3/(r1^3))*(vv1+vv1');
%v2 = cross(s2,n2);
%vv2 = n2'*v2;
%As2 = (3/(r2^3))*(vv2+vv2');
%A = Ap1+Ap2; % +As1+As2;
%Asq = sum(sum(A.^2));
%Asq; % output
% syms x y z p c r1 r2
% c1 = [c,0,0];
% c2 = -1*c1;
% p1 = [0,p,0];
% p2 = -1*p1;
% n1 = [x-c1(1),y,z]/r1;
% n2 = [x-c2(1),y,z]/r2;
% mew = eye(3);
% Ap1 = (p1'*n1+n1'*p1-(mew-n1'*n1)*dot(n1,p1)) *(3/(2*(r1)^2))
% Ap2 = (p2'*n2+n2'*p2-(mew-n2'*n2)*dot(n2,p2)) *(3/(2*(r2)^2))

