function problem1
% löst -u'' + u = 1
a = 0;
b = 1;
ul = 0;
ur = 0;

n = 100;
h = 1/(n-1);
x = linspace(a, b, n);

m = 1/6 * [2 1
           1 2];
k = [ 1 -1
     -1  1];

B = zeros(n-1, 2);

parfor i=1:n-1
    B(i, :) = [i, i+1];
end

M = createMatrix(m, B, h);
K = createMatrix(k, B, 1/h);

A = M+K;

b = h/2 * [1; 2 * ones(n-2, 1); 1];

[A, b] = neumannBoundary(A, b, ul, ur);

u = A\b;

factor = -1/(1+exp(1));
y = factor*exp(x)+(-factor-1)*exp(-x)+1;

plot(x, u, x, y);
end