function fem
a = 0;
b = 5;
n = 6;
x = linspace(a,b,n);
h = (b-a)/(n-1);

f = @(x) ones(1,length(x));

for i=1:n
    phi(i) = LinearBasisFunction(x(i) - h, x(i), x(i) + h);
end

M = zeros(n,n);
K = zeros(n,n);
bvector = zeros(n,1);

for i=1:n
    fv = @(x) f(x) .* phi(i).evaluate(x);
    
    bvector(i) = quadgk(fv, a, b);
    
    for j = 1:n
        product = @(x) phi(i).evaluate(x) .* phi(j).evaluate(x);
        diff = @(x) phi(i).evaluateDerivate(x) .* phi(j).evaluateDerivate(x);
        
        M(i,j) = quadgk(product, a, b);
        K(i,j) = quadgk(diff, a, b);
    end
end

K
M
bvector

u=bvector\(K+M);

f=-1/(1+exp(5));
y = f*exp(x)+(-f-1)*exp(-x)+1;

plot(x, u, x, y);

end