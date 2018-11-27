N = 8;

% C_ijk = A_pi * B_pjk
A = rand(N,N);
B = rand(N,N,N);
C = zeros(size(B));
for j = 1:N
    C(:,:,j) = A*B(:,:,j);
end

% C_ijk = A_qj * B_iqk
A = rand(N,N);
B = rand(N,N,N);
B = permute(B,[]);
C = zeros(size(B));
for j = 1:N
    C(:,:,j) = A*B(:,:,j);
end

% C_ijk = A_rk * B_ijr
A = rand(N,N);
B = rand(N,N,N);
B = permute(B,[]);
C = zeros(size(B));
for j = 1:N
    C(:,:,j) = A*B(:,:,j);
end