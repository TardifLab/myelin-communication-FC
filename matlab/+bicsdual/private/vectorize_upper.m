function v = vectorize_upper(A)
% Upper triangle (i<j).
n = size(A,1);
mask = triu(true(n),1);
v = A(mask);
end
