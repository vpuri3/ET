p = 2; q = 5; r = p*q;
I = eye(r);
c = 1:q;
S = [];
for c = 1:q
    S = [S; I(c:q:r,:)];
end
S