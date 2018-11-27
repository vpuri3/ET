function v = mat2vec(A)
s = size(A);
d = length(s); % number of dimensions.
n1 = s(1); n2 = s(2);
v = zeros(n1*n2,1);
for d = 1:n2 %traversing y
    for c = 1:n1 %traversing x
        k = c+n2*(d-1);
        v(k) = A(c,d);
    end
end

end
