function fem
a = 0;
b = 3;
ul = 0;
ur = 3;
n = 31;
x = linspace(a,b,n);
h = (b-a)/(n-1);

f = @(x) ones(1,length(x));

phi(1) = LinearBasisFunction(-1,0,1);
phi(2) = LinearBasisFunction(0,1,2);

for i=1:n
    phiGlob(i) = LinearBasisFunction(x(i)-h,x(i),x(i)+h);
end

l = length(phi);
m = zeros(l);
k = zeros(l);

for i=1:l    
    for j=1:l
        product = @(x) phi(i).evaluate(x) .* phi(j).evaluate(x);
        diffProduct = @(x) phi(i).evaluateDerivate(x) .* phi(j).evaluateDerivate(x);
        
        m(i,j) = quadgk(product, 0, 1)*h;
        k(i,j) = quadgk(diffProduct, 0, 1)/h;
    end
end

bvector = zeros(n,1);
M = zeros(n,n);
K = zeros(n,n);

for i=1:n    
    bvector(i) = quadgk(@(x) f(x) .* phiGlob(i).evaluate(x), a, b);
        
    parfor j = 1:n
        if i-j == 1
            M(i,j) = m(1,2);
            K(i,j) = k(1,2);
        elseif i-j == 0
            if i == 1 || i == n
                M(i,j) = m(1,1); 
                K(i,j) = k(1,1);
            else
                M(i,j) = m(1,1) + m(2,2);
                K(i,j) = k(1,1) + k(2,2);
            end
        elseif i-j == -1
            M(i,j) = m(2,1);
            K(i,j) = k(2,1);
        end
    end
end

A = K+M;

bvector = bvector-A(:,1)*ul-A(:,n)*ur;
bvector(1) = ul;
bvector(n) = ur;

A(1,2:n) = zeros(1,n-1);
A(:,1) = zeros(n,1);
A(:,n) = zeros(n,1);
A(n,2:n) = zeros(1,n-1);
A(1,1) = 1;
A(n,n) = 1;

u=A\bvector;

f = -1/(1+exp(3));
y = f*exp(x)+(-f-1)*exp(-x)+1;

plot(x, u, x, y);

end