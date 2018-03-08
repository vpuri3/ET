% Vx is matrix of eigenvectors
for c = 1:n
    s = Ax*V(:,c);
    s = s ./ V(:,c)
    Dx(c,c)
end