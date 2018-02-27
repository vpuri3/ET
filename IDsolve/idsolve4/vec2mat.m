function A = vec2mat(v,m,n)
d = length(v);
assert(m*n == d);
A = zeros(m,n);
for c = 1:n
    k = (c-1)*m + 1;
    A(:,c) = v(k:k+m-1);
end
